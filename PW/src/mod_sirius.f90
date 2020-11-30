MODULE mod_sirius
USE input_parameters, ONLY : use_sirius, sirius_cfg, use_sirius_scf
USE sirius
USE funct
IMPLICIT NONE

! use SIRIUS to solve KS equations
LOGICAL :: use_sirius_ks_solver               = .TRUE.
! use SIRIUS to generate density
LOGICAL :: use_sirius_density                 = .TRUE.
! use SIRIUS to generate effective potential; WARNING: currently must be always set to .false.
LOGICAL :: use_sirius_potential               = .FALSE.
! use SIRIUS to generate density matrix ('bec' thing in QE) WARNING: currently must be set to the value of use_sirius_density
LOGICAL :: use_sirius_density_matrix          = .TRUE.
! use SIRIUS to compute local part of pseudopotential
LOGICAL :: use_sirius_vloc                    = .TRUE.
! use SIRIUS to compute core charge density
LOGICAL :: use_sirius_rho_core                = .TRUE.
! use SIRIUS to compute plane-wave coefficients of atomic charge density
LOGICAL :: use_sirius_rho_atomic              = .TRUE.
! use SIRIUS to compute forces
LOGICAL :: use_sirius_forces                  = .TRUE.
! use SIRIUS to compute stress tensor
LOGICAL :: use_sirius_stress                  = .TRUE.

! inverse of the reciprocal lattice vectors matrix
REAL(8) bg_inv(3,3)
! total number of k-points
INTEGER num_kpoints
REAL(8), ALLOCATABLE :: kpoints(:,:)
REAL(8), ALLOCATABLE :: wkpoints(:)
REAL(8), ALLOCATABLE :: beta_ri_tab(:,:,:)
REAL(8), ALLOCATABLE :: aug_ri_tab(:,:,:,:)

TYPE atom_type_t
  ! atom label
  CHARACTER(len=100, kind=C_CHAR) :: label
  ! nh(iat) in the QE notation
  INTEGER                 :: num_beta_projectors
  ! lmax for beta-projectors
  INTEGER                 :: lmax
  ! plane-wave coefficients of Q-operator
  COMPLEX(8), ALLOCATABLE :: qpw(:, :)
END TYPE atom_type_t

TYPE(atom_type_t), ALLOCATABLE :: atom_type(:)

TYPE(C_PTR) :: sctx = C_NULL_PTR
TYPE(C_PTR) :: gs_handler = C_NULL_PTR
TYPE(C_PTR) :: ks_handler = C_NULL_PTR

CONTAINS

! A callback function to compute band occupancies is SIRIUS using QE
! TODO: merge to a single function with argments (energies in, occupancies out)
SUBROUTINE calc_band_occupancies() BIND(C)
USE ISO_C_BINDING
IMPLICIT NONE
!
CALL get_band_energies_from_sirius()
CALL weights()
CALL put_band_occupancies_to_sirius()
!
END SUBROUTINE calc_band_occupancies


! A callback function to compute radial integrals of Vloc(r)
SUBROUTINE calc_vloc_radial_integrals(iat, nq, q, vloc_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE esm
USE Coul_cut_2D
USE uspp_param, ONLY : upf
IMPLICIT NONE
!
INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
REAL(KIND=c_double), INTENT(IN)  :: q(nq)
REAL(KIND=c_double), INTENT(OUT) :: vloc_ri(nq)
!
REAL(DP) :: aux(rgrid(iat)%mesh), aux1(rgrid(iat)%mesh), r, q2
INTEGER :: iq, ir, nr
!
REAL(DP), EXTERNAL :: qe_erf
! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
nr = msh(iat)
! q-independent radial function
DO ir = 1, nr
  r = rgrid(iat)%r(ir)
  aux1(ir) = r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * qe_erf(r)
END DO
! loop over q-points
DO iq = 1, nq
  IF (q(iq) < eps8) THEN ! q=0 case
    !
    ! first the G=0 term
    !
    IF ((do_comp_esm .AND. (esm_bc .NE. 'pbc')) .OR. do_cutoff_2D) THEN
      !
      ! ... temporarily redefine term for ESM calculation
      !
      DO ir = 1, nr
        r = rgrid(iat)%r(ir)
        aux(ir) = r * (r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * qe_erf(r))
      END DO
      IF (do_cutoff_2D .AND. rgrid(iat)%r(nr) > lz) THEN
        CALL errore('vloc_of_g','2D cutoff is smaller than pseudo cutoff radius: &
          & increase interlayer distance (or see Modules/read_pseudo.f90)',1)
      END IF
    ELSE
      DO ir = 1, nr
        r = rgrid(iat)%r(ir)
        aux(ir) = r * (r * upf(iat)%vloc(ir) + upf(iat)%zp * e2)
      END DO
    END IF
    CALL simpson(nr, aux, rgrid(iat)%rab(1), vloc_ri(iq))
  ELSE ! q > 0 case
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = aux1(ir) * SIN(q(iq) * r) / q(iq)
    END DO
    CALL simpson(nr, aux, rgrid(iat)%rab(1), vloc_ri(iq))
    IF ((.NOT.do_comp_esm) .OR. (esm_bc .EQ. 'pbc')) THEN
      !
      !   here we re-add the analytic fourier transform of the erf function
      !
      IF (.NOT. do_cutoff_2D) THEN
        q2 = q(iq) * q(iq)
        vloc_ri(iq) = vloc_ri(iq) - upf(iat)%zp * e2 * EXP(-q2 * 0.25d0) / q2
      END IF
    END IF
  END IF
  ! convert to Ha
  vloc_ri(iq) = vloc_ri(iq) / 2.0
END DO

END SUBROUTINE calc_vloc_radial_integrals


! A callback function to compute radial integrals of Vloc(r) with
! derivatives of Bessel functions
SUBROUTINE calc_vloc_dj_radial_integrals(iat, nq, q, vloc_dj_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE esm
USE Coul_cut_2D
USE uspp_param, ONLY : upf
IMPLICIT NONE
!
INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
REAL(KIND=c_double), INTENT(IN)  :: q(nq)
REAL(KIND=c_double), INTENT(OUT) :: vloc_dj_ri(nq)
!
REAL(DP) :: aux(rgrid(iat)%mesh), aux1(rgrid(iat)%mesh), r, q2
INTEGER :: iq, ir, nr
!
REAL(DP), EXTERNAL :: qe_erf
! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
nr = msh(iat)
! q-independent radial function
DO ir = 1, nr
  r = rgrid(iat)%r(ir)
  aux1(ir) = r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * qe_erf(r)
END DO
! loop over q-points
DO iq = 1, nq
  IF (q(iq) < eps8) THEN ! q=0 case
    vloc_dj_ri(iq) = 0.d0
  ELSE ! q > 0 case
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = aux1(ir) * (SIN(q(iq) * r) / q(iq)**2 - r * COS(q(iq) * r) / q(iq))
    END DO
    CALL simpson(nr, aux, rgrid(iat)%rab(1), vloc_dj_ri(iq))
    vloc_dj_ri(iq) = vloc_dj_ri(iq) / q(iq)
    IF ((.NOT.do_comp_esm) .OR. (esm_bc .EQ. 'pbc')) THEN
      IF (.NOT. do_cutoff_2D) THEN
        q2 = q(iq) * q(iq)
        vloc_dj_ri(iq) = vloc_dj_ri(iq) - upf(iat)%zp * e2 * &
          &EXP(-q2 * 0.25d0) * (q2 + 4) / 2 / q2 / q2
      END IF
    END IF
  END IF
  ! convert to Ha
  vloc_dj_ri(iq) = vloc_dj_ri(iq) / 2.0
END DO

END SUBROUTINE calc_vloc_dj_radial_integrals


! A callback function to compute radial integrals of rho_core(r)
SUBROUTINE calc_rhoc_radial_integrals(iat, nq, q, rhoc_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE uspp_param, ONLY : upf
IMPLICIT NONE
!
INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
REAL(KIND=c_double), INTENT(IN)  :: q(nq)
REAL(KIND=c_double), INTENT(OUT) :: rhoc_ri(nq)
!
REAL(DP) :: aux(rgrid(iat)%mesh), r
INTEGER :: iq, ir, nr
! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
nr = msh(iat)
! loop over q-points
DO iq = 1, nq
  IF (q(iq) < eps8) THEN ! q=0 case
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = r**2 * upf(iat)%rho_atc(ir)
    ENDDO
  ELSE
    CALL sph_bes(nr, rgrid(iat)%r(1), q(iq), 0, aux)
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = r**2 * upf(iat)%rho_atc(ir) * aux(ir)
    ENDDO
  END IF
  CALL simpson(nr, aux, rgrid(iat)%rab(1), rhoc_ri(iq))
END DO

END SUBROUTINE calc_rhoc_radial_integrals


! A callbacl function to compute radial integrals or rho_core(r) with
! the derivatives of Bessel functions
SUBROUTINE calc_rhoc_dj_radial_integrals(iat, nq, q, rhoc_dj_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE uspp_param, ONLY : upf
IMPLICIT NONE
!
INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
REAL(KIND=c_double), INTENT(IN)  :: q(nq)
REAL(KIND=c_double), INTENT(OUT) :: rhoc_dj_ri(nq)
!
REAL(DP) :: aux(rgrid(iat)%mesh), r
INTEGER :: iq, ir, nr
! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
nr = msh(iat)
! loop over q-points
DO iq = 1, nq
  IF (q(iq) < eps8) THEN ! q=0 case
    rhoc_dj_ri(iq) = 0.d0
  ELSE
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = r * upf(iat)%rho_atc(ir) * &
        &(r * COS(q(iq) * r) / q(iq) - SIN(q(iq) * r) / q(iq)**2)
    ENDDO
    CALL simpson(nr, aux, rgrid(iat)%rab(1), rhoc_dj_ri(iq))
  END IF
END DO

END SUBROUTINE calc_rhoc_dj_radial_integrals


SUBROUTINE calc_beta_radial_integrals(iat, q, beta_ri, ld) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, tab
USE uspp_param,   ONLY : upf
IMPLICIT NONE
!
INTEGER(kind=c_int), INTENT(in), VALUE :: iat
REAL(kind=c_double), INTENT(in), VALUE :: q
INTEGER(kind=c_int), INTENT(in), VALUE :: ld
REAL(kind=c_double), INTENT(out) :: beta_ri(ld)
!
REAL(8) :: px, ux, vx, wx
INTEGER :: i0, i1, i2, i3, ib
!
IF (ld.LT.upf(iat)%nbeta) THEN
  WRITE(*,*)'not enough space to store all beta projectors, ld=',ld,' nbeta=',upf(iat)%nbeta
  STOP
ENDIF
IF (.NOT.ALLOCATED(tab)) THEN
  WRITE(*,*)'tab array is not allocated'
  STOP
ENDIF

px = q / dq - INT(q / dq)
ux = 1.d0 - px
vx = 2.d0 - px
wx = 3.d0 - px
i0 = INT(q / dq) + 1
i1 = i0 + 1
i2 = i0 + 2
i3 = i0 + 3
DO ib = 1, upf(iat)%nbeta
  beta_ri(ib) = beta_ri_tab(i0, ib, iat) * ux * vx * wx / 6.d0 + &
                beta_ri_tab(i1, ib, iat) * px * vx * wx / 2.d0 - &
                beta_ri_tab(i2, ib, iat) * px * ux * wx / 2.d0 + &
                beta_ri_tab(i3, ib, iat) * px * ux * vx / 6.d0
ENDDO

END SUBROUTINE calc_beta_radial_integrals


SUBROUTINE calc_beta_dj_radial_integrals(iat, q, beta_ri, ld) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, tab
USE uspp_param,   ONLY : upf
IMPLICIT NONE
!
INTEGER(kind=c_int), INTENT(in), VALUE :: iat
REAL(kind=c_double), INTENT(in), VALUE :: q
INTEGER(kind=c_int), INTENT(in), VALUE :: ld
REAL(kind=c_double), INTENT(out) :: beta_ri(ld)
!
REAL(8) :: px, ux, vx, wx
INTEGER :: i0, i1, i2, i3, ib
!
IF (ld.LT.upf(iat)%nbeta) THEN
  WRITE(*,*)'not enough space to store all beta projectors, ld=',ld,' nbeta=',upf(iat)%nbeta
  STOP
ENDIF
IF (.NOT.ALLOCATED(tab)) THEN
  WRITE(*,*)'tab array is not allocated'
  STOP
ENDIF

px = q / dq - INT(q / dq)
ux = 1.d0 - px
vx = 2.d0 - px
wx = 3.d0 - px
i0 = INT(q / dq) + 1
i1 = i0 + 1
i2 = i0 + 2
i3 = i0 + 3
DO ib = 1, upf(iat)%nbeta
  beta_ri(ib) = beta_ri_tab(i0, ib, iat) * (-vx*wx-ux*wx-ux*vx)/6.d0 + &
                beta_ri_tab(i1, ib, iat) * (+vx*wx-px*wx-px*vx)/2.d0 - &
                beta_ri_tab(i2, ib, iat) * (+ux*wx-px*wx-px*ux)/2.d0 + &
                beta_ri_tab(i3, ib, iat) * (+ux*vx-px*vx-px*ux)/6.d0
  beta_ri(ib) = beta_ri(ib) / dq
ENDDO

END SUBROUTINE calc_beta_dj_radial_integrals


SUBROUTINE calc_aug_radial_integrals(iat, q, aug_ri, ld1, ld2) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, qrad
USE uspp_param,   ONLY : upf
IMPLICIT NONE
!
INTEGER(kind=c_int), INTENT(in), VALUE :: iat
REAL(kind=c_double), INTENT(in), VALUE :: q
INTEGER(kind=c_int), INTENT(in), VALUE :: ld1
INTEGER(kind=c_int), INTENT(in), VALUE :: ld2
REAL(kind=c_double), INTENT(out) :: aug_ri(ld1, ld2)
!
REAL(8) :: px, ux, vx, wx
INTEGER :: i0, i1, i2, i3, l, nb, mb, ijv
!
IF (upf(iat)%tvanp) THEN
  px = q / dq - INT(q / dq)
  ux = 1.d0 - px
  vx = 2.d0 - px
  wx = 3.d0 - px
  i0 = INT(q / dq) + 1
  i1 = i0 + 1
  i2 = i0 + 2
  i3 = i0 + 3
  DO l = 1, 2 * atom_type(iat)%lmax + 1
    DO nb = 1, upf(iat)%nbeta
      DO mb = nb, upf(iat)%nbeta
        ijv = mb * (mb-1) / 2 + nb
        aug_ri(ijv, l) = aug_ri_tab(i0, ijv, l, iat) * ux * vx * wx / 6.d0 + &
                         aug_ri_tab(i1, ijv, l, iat) * px * vx * wx / 2.d0 - &
                         aug_ri_tab(i2, ijv, l, iat) * px * ux * wx / 2.d0 + &
                         aug_ri_tab(i3, ijv, l, iat) * px * ux * vx / 6.d0
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE calc_aug_radial_integrals


SUBROUTINE calc_aug_dj_radial_integrals(iat, q, aug_ri, ld1, ld2) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, qrad
USE uspp_param,   ONLY : upf
IMPLICIT NONE
!
INTEGER(kind=c_int), INTENT(in), VALUE :: iat
REAL(kind=c_double), INTENT(in), VALUE :: q
INTEGER(kind=c_int), INTENT(in), VALUE :: ld1
INTEGER(kind=c_int), INTENT(in), VALUE :: ld2
REAL(kind=c_double), INTENT(out) :: aug_ri(ld1, ld2)
!
REAL(8) :: px, ux, vx, wx
INTEGER :: i0, i1, i2, i3, l, nb, mb, ijv
!
IF (upf(iat)%tvanp) THEN
  px = q / dq - INT(q / dq)
  ux = 1.d0 - px
  vx = 2.d0 - px
  wx = 3.d0 - px
  i0 = INT(q / dq) + 1
  i1 = i0 + 1
  i2 = i0 + 2
  i3 = i0 + 3
  DO l = 1, 2 * atom_type(iat)%lmax + 1
    DO nb = 1, upf(iat)%nbeta
      DO mb = nb, upf(iat)%nbeta
        ijv = mb * (mb-1) / 2 + nb
        aug_ri(ijv, l) = - aug_ri_tab(i0, ijv, l, iat) * (ux*vx + vx*wx + ux*wx) / 6.d0 &
                         + aug_ri_tab(i1, ijv, l, iat) * (wx*vx - px*wx - px*vx) * 0.5d0 &
                         - aug_ri_tab(i2, ijv, l, iat) * (wx*ux - px*wx - px*ux) * 0.5d0 &
                         + aug_ri_tab(i3, ijv, l, iat) * (ux*vx - px*ux - px*vx) / 6.d0
        aug_ri(ijv, l) = aug_ri(ijv, l) / dq
      ENDDO
    ENDDO
  ENDDO
ENDIF

END SUBROUTINE calc_aug_dj_radial_integrals


SUBROUTINE setup_sirius()
USE cell_base, ONLY : alat, at, bg
USE funct, ONLY : get_iexch, get_icorr, get_inlc, get_meta, get_igcc, get_igcx
USE ions_base, ONLY : tau, nsp, atm, zv, amass, ityp, nat
USE uspp_param, ONLY : upf, nhm, nh
USE atom, ONLY : rgrid, msh
USE fft_base, ONLY :  dfftp
USE klist, ONLY : nks, xk, nkstot, wk
USE gvect, ONLY : ngm_g, ecutrho, ngm, mill
USE gvecw, ONLY : ecutwfc
USE control_flags, ONLY : gamma_only, diago_full_acc, mixing_beta, nmix
USE mp_pools, ONLY : inter_pool_comm, npool
USE mp_images,        ONLY : nproc_image, intra_image_comm
USE mp, ONLY : mp_sum, mp_bcast
USE wvfct, ONLY : nbnd
USE parallel_include
USE sirius
USE input_parameters, ONLY : sirius_cfg, diago_david_ndim
USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
USE lsda_mod, ONLY : lsda, nspin, starting_magnetization
USE cell_base, ONLY : omega
USE symm_base, ONLY : nosym
USE spin_orb,  ONLY : lspinorb
USE ldaU, ONLY : lda_plus_U, Hubbard_J, Hubbard_U, Hubbard_alpha, &
     & Hubbard_beta, is_Hubbard, lda_plus_u_kind, Hubbard_J0, U_projection, Hubbard_l
USE esm,       ONLY : esm_local, esm_bc, do_comp_esm
USE control_flags, ONLY : iverbosity
USE Coul_cut_2D,          ONLY : do_cutoff_2D
IMPLICIT NONE
!
INTEGER :: dims(3), i, ia, iat, rank, ierr, ijv, li, lj, mb, nb, j, l,&
     ilast, ir, num_gvec, num_ranks_k, vt(3), iwf, num_kp, nmagd
REAL(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp
REAL(8), ALLOCATABLE :: dion(:, :), qij(:,:,:), vloc(:), wk_tmp(:), xk_tmp(:,:)
INTEGER, ALLOCATABLE :: nk_loc(:)
INTEGER :: ih, jh, ijh, lmax_beta
CHARACTER(LEN=1024) :: conf_str
INTEGER, EXTERNAL :: set_hubbard_l,set_hubbard_n
REAL(8), EXTERNAL :: hubbard_occ


! TODO: check if this is necessary now, them radial integrals of vloc are
! computed by QE
IF (do_cutoff_2D) THEN
  use_sirius_vloc = .FALSE.
ENDIF

ALLOCATE(atom_type(nsp))

DO iat = 1, nsp
  atom_type(iat)%label = TRIM(ADJUSTL(atm(iat)))

  lmax_beta = 0
  DO i = 1, upf(iat)%nbeta
    lmax_beta = MAX(lmax_beta, upf(iat)%lll(i))
  ENDDO
  atom_type(iat)%lmax = lmax_beta
ENDDO

! create context of simulation
CALL sirius_create_context(intra_image_comm, sctx)
! create initial configuration dictionary in JSON
WRITE(conf_str, 10)diago_david_ndim, mixing_beta, nmix+1
10 FORMAT('{"parameters" : {"electronic_structure_method" : "pseudopotential"},&
            &"iterative_solver" : {"residual_tolerance" : 1e-6,"subspace_size" : ',I4,'}, &
            &"mixer" : {"beta" : ',F12.6,',"max_history" : ',I4,', "use_hartree" : true},&
            &"settings" : {"itsol_tol_scale" : [0.1, 1.0]}}')
! set initial parameters
CALL sirius_import_parameters(sctx, conf_str)
! set default verbosity
CALL sirius_set_parameters(sctx, verbosity=MIN(1, iverbosity))
! import config file
CALL sirius_import_parameters(sctx, TRIM(ADJUSTL(sirius_cfg)))

! derive the number of magnetic dimensions
nmagd = 0
IF (nspin.EQ.2) THEN
  nmagd = 1
ENDIF
IF (lspinorb.OR.noncolin) THEN
  nmagd = 3
ENDIF
! copy dimensions of dense FFT grid
dims(1) = dfftp%nr1
dims(2) = dfftp%nr2
dims(3) = dfftp%nr3
! set basic parameters
! set |G| cutoff of the dense FFT grid: convert from G^2/2 Rydbergs to |G| in [a.u.^-1]
! set |G+k| cutoff for the wave-functions: onvert from |G+k|^2/2 Rydbergs to |G+k| in [a.u.^-1]
! use symmetrization either on SIRIUS or QE side
CALL sirius_set_parameters(sctx, num_bands=nbnd, num_mag_dims=nmagd, gamma_point=gamma_only,&
  &use_symmetry=use_sirius_scf.AND..NOT.nosym, so_correction=lspinorb,&
  &pw_cutoff=SQRT(ecutrho), gk_cutoff=SQRT(ecutwfc),&
  &hubbard_correction=lda_plus_U, hubbard_correction_kind=lda_plus_u_kind,&
  &hubbard_orbitals=TRIM(ADJUSTL(U_projection)),fft_grid_size=dims)

! check if this is requred now then radial integrals of Vloc are computed by QE
IF (do_comp_esm) THEN
  CALL sirius_set_parameters(sctx, esm_bc=esm_bc)
ENDIF

num_ranks_k = nproc_image / npool
i = SQRT(DBLE(num_ranks_k) + 1d-10)
IF (i * i .NE. num_ranks_k) THEN
  dims(1) = 1
  dims(2) = num_ranks_k
ELSE
  dims(1) = i
  dims(2) = i
ENDIF
CALL sirius_set_mpi_grid_dims(sctx, 2, dims(1))

IF (diago_full_acc) THEN
  CALL sirius_set_parameters(sctx, iter_solver_tol_empty=0.d0)
ELSE
  CALL sirius_set_parameters(sctx, iter_solver_tol_empty=1d-5)
ENDIF

CALL sirius_set_callback_function(sctx, "beta_ri", C_FUNLOC(calc_beta_radial_integrals))
CALL sirius_set_callback_function(sctx, "beta_ri_djl", C_FUNLOC(calc_beta_dj_radial_integrals))
CALL sirius_set_callback_function(sctx, "aug_ri", C_FUNLOC(calc_aug_radial_integrals))
CALL sirius_set_callback_function(sctx, "aug_ri_djl", C_FUNLOC(calc_aug_dj_radial_integrals))
CALL sirius_set_callback_function(sctx, "vloc_ri", C_FUNLOC(calc_vloc_radial_integrals))
CALL sirius_set_callback_function(sctx, "vloc_ri_djl", C_FUNLOC(calc_vloc_dj_radial_integrals))
CALL sirius_set_callback_function(sctx, "rhoc_ri", C_FUNLOC(calc_rhoc_radial_integrals))
CALL sirius_set_callback_function(sctx, "rhoc_ri_djl", C_FUNLOC(calc_rhoc_dj_radial_integrals))
CALL sirius_set_callback_function(sctx, "band_occ", C_FUNLOC(calc_band_occupancies))

!call sirius_set_parameters(sctx, min_occupancy=0.01d0)

! set lattice vectors of the unit cell (length is in [a.u.])
a1(:) = at(:, 1) * alat
a2(:) = at(:, 2) * alat
a3(:) = at(:, 3) * alat
CALL sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))

vlat(:, 1) = a1(:)
vlat(:, 2) = a2(:)
vlat(:, 3) = a3(:)
! get the inverse of Bravais lattice vectors
CALL invert_mtrx(vlat, vlat_inv)
! get the inverse of reciprocal lattice vectors
CALL invert_mtrx(bg, bg_inv)
! get MPI rank associated with the distribution of k-points
CALL mpi_comm_rank(inter_pool_comm, rank, ierr)

! initialize atom types
DO iat = 1, nsp

  ! add new atom type
   CALL sirius_add_atom_type(sctx, atom_type(iat)%label, &
        & zn=NINT(zv(iat)+0.001d0), &
        & mass=amass(iat), &
        & spin_orbit=upf(iat)%has_so)


  ! set radial grid
  CALL sirius_set_atom_type_radial_grid(sctx, atom_type(iat)%label, upf(iat)%mesh, upf(iat)%r)

  ! set beta-projectors
  DO i = 1, upf(iat)%nbeta
    l = upf(iat)%lll(i);
    IF (upf(iat)%has_so) THEN
      IF (upf(iat)%jjj(i) .LE. upf(iat)%lll(i)) THEN
        l = - upf(iat)%lll(i)
      ENDIF
    ENDIF
    CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "beta", &
         & upf(iat)%beta(1:upf(iat)%kbeta(i), i), upf(iat)%kbeta(i), l=l)
  ENDDO

  ! set the atomic radial functions
  DO iwf = 1, upf(iat)%nwfc
    l = upf(iat)%lchi(iwf)
    IF (upf(iat)%has_so) THEN
      IF (upf(iat)%jchi(iwf) < l) THEN
        l = -l
      ENDIF
    ENDIF
    IF (ALLOCATED(upf(iat)%nchi)) THEN
      i = upf(iat)%nchi(iwf)
    ELSE
      i = -1
    ENDIF
    CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "ps_atomic_wf", &
         & upf(iat)%chi(1:msh(iat), iwf), msh(iat), l=l, occ=upf(iat)%oc(iwf), n=i)
  ENDDO

  IF (is_hubbard(iat)) THEN
     CALL sirius_set_atom_type_hubbard(sctx, &
          & atom_type(iat)%label, &
          & Hubbard_l(iat), &
          & set_hubbard_n(upf(iat)%psd), &
          & hubbard_occ ( upf(iat)%psd ), &
          & Hubbard_U(iat), &
          & Hubbard_J(1,iat), &
          & Hubbard_alpha(iat), &
          & Hubbard_beta(iat), &
          & Hubbard_J0(iat))
  ENDIF

  ALLOCATE(dion(upf(iat)%nbeta, upf(iat)%nbeta))
  ! convert to hartree
  DO i = 1, upf(iat)%nbeta
    DO j = 1, upf(iat)%nbeta
      dion(i, j) = upf(iat)%dion(i, j) / 2.d0
    END DO
  END DO
  ! sed d^{ion}_{i,j}
  CALL sirius_set_atom_type_dion(sctx, atom_type(iat)%label, upf(iat)%nbeta, dion(1, 1))
  DEALLOCATE(dion)

  ! get lmax_beta for this atom type
  lmax_beta = -1
  DO i = 1, upf(iat)%nbeta
    lmax_beta = MAX(lmax_beta, upf(iat)%lll(i))
  ENDDO

  ! set radial function of augmentation charge
  IF (upf(iat)%tvanp) THEN
    !do l = 0, upf(iat)%nqlc - 1
    DO l = 0, 2 * lmax_beta
      DO i = 0, upf(iat)%nbeta - 1
        DO j = i, upf(iat)%nbeta - 1
          ijv = j * (j + 1) / 2 + i + 1
          CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "q_aug",&
                                                   &upf(iat)%qfuncl(1:upf(iat)%kkbeta, ijv, l), upf(iat)%kkbeta,&
                                                   &l=l, idxrf1=i, idxrf2=j)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF (upf(iat)%tpawp) THEN
    DO i = 1, upf(iat)%nbeta
      CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "ae_paw_wf",&
                                               &upf(iat)%aewfc(1:upf(iat)%paw%iraug,i), upf(iat)%paw%iraug)
      CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "ps_paw_wf",&
                                               &upf(iat)%pswfc(1:upf(iat)%paw%iraug,i), upf(iat)%paw%iraug)
    ENDDO
    CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "ae_paw_core",&
                                             &upf(iat)%paw%ae_rho_atc, upf(iat)%mesh)

    CALL sirius_set_atom_type_paw(sctx, atom_type(iat)%label, upf(iat)%paw%core_energy / 2,&
                                 &upf(iat)%paw%oc, upf(iat)%nbeta)
  ENDIF

  ! set non-linear core correction
  IF (.TRUE.) THEN !use_sirius_rho_core) THEN
    ALLOCATE(vloc(upf(iat)%mesh))
    vloc = 0.d0
    IF (ALLOCATED(upf(iat)%rho_atc)) THEN
      DO i = 1, msh(iat)
        vloc(i) = upf(iat)%rho_atc(i)
      ENDDO
    ENDIF
    CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "ps_rho_core",&
                                             &vloc, upf(iat)%mesh)
    DEALLOCATE(vloc)
  ENDIF

  ! set total charge density of a free atom (to compute initial rho(r))
  CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "ps_rho_total",&
                                           &upf(iat)%rho_at, upf(iat)%mesh)

  ! the hack is done in Modules/readpp.f90
  IF (.TRUE.) THEN !use_sirius_vloc) THEN
    ALLOCATE(vloc(upf(iat)%mesh))
    DO i = 1, msh(iat)
      vloc(i) = upf(iat)%vloc(i)
    ENDDO
    ! convert to Hartree
    vloc = vloc / 2.d0
    ! add a correct tail
    DO i = msh(iat) + 1, upf(iat)%mesh
      vloc(i) = -zv(iat) / upf(iat)%r(i)
    ENDDO
    ! set local part of pseudo-potential
    CALL sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, "vloc", vloc, upf(iat)%mesh)
    DEALLOCATE(vloc)
  ENDIF
ENDDO

! add atoms to the unit cell
! WARNING: sirius accepts only fractional coordinates;
!          if QE stores coordinates in a different way, the conversion must be made here
DO ia = 1, nat
  iat = ityp(ia)
  ! Cartesian coordinates
  v1(:) = tau(:, ia) * alat
  ! fractional coordinates
  v1(:) = MATMUL(vlat_inv, v1)
  ! reduce coordinates to [0, 1) interval
  !call sirius_reduce_coordinates(v1(1), v2(1), vt(1))
  v2 = v1
  IF (noncolin) THEN
    v1(1) = zv(iat) * starting_magnetization(iat) * SIN(angle1(iat)) * COS(angle2(iat))
    v1(2) = zv(iat) * starting_magnetization(iat) * SIN(angle1(iat)) * SIN(angle2(iat))
    v1(3) = zv(iat) * starting_magnetization(iat) * COS(angle1(iat))
  ELSE
    v1 = 0
    v1(3) = zv(iat) * starting_magnetization(iat)
  ENDIF
  CALL sirius_add_atom(sctx, atom_type(iat)%label, v2, v1)
ENDDO

CALL put_xc_functional_to_sirius()

write(*,*)''
write(*,*)'=========================================='
write(*,*)'* initializing SIRIUS simulation context *'
write(*,*)'=========================================='

! initialize global variables/indices/arrays/etc. of the simulation
CALL sirius_initialize_context(sctx)

DO iat = 1, nsp
  CALL sirius_get_num_beta_projectors(sctx, atom_type(iat)%label, atom_type(iat)%num_beta_projectors)
ENDDO

!! get number of g-vectors of the dense fft grid
!call sirius_get_num_gvec(num_gvec)
!
!! TODO: number of G-vectors can be different; adapt the code wo work in this situation
!if (.not.((num_gvec .eq. ngm_g) .or. (num_gvec * 2 - 1 .eq. ngm_g))) then
!  write(*,*)"wrong number of g-vectors"
!  write(*,*)"num_gvec=",num_gvec
!  write(*,*)"ngm_g=",ngm_g
!endif

!call sirius_get_fft_grid_size(dims(1)) ! TODO: size of FFT box is not very relevant and in principle can be slightly different
!if (dims(1).ne.dfftp%nr1.or.dims(2).ne.dfftp%nr2.or.dims(3).ne.dfftp%nr3) then
!  write(*,*)"wrong fft grid dimensions"
!  write(*,*)"qe: ", dfftp%nr1,  dfftp%nr2,  dfftp%nr3
!  write(*,*)"sirius: ", dims
!  stop 111
!endif

!if (.true.) then
!  do i = 1, num_kpoints
!    write(*,*)'ik=',i,' kpoint=',matmul(bg_inv,kpoints(:,i))
!  enddo
!endif

!allocate(wk_tmp(nkstot))
!allocate(xk_tmp(3, nkstot))
!! weights of k-points in SIRIUS must sum to one
!do i = 1, nkstot
!  if (nspin.eq.1) then
!    wk_tmp(i) = wk(i) / 2.d0
!  else
!    wk_tmp(i) = wk(i)
!  endif
!  xk_tmp(:,i) = xk(:,i)
!end do
!
!call mpi_bcast(wk_tmp(1),        nkstot, mpi_double, 0, inter_pool_comm, ierr)
!call mpi_bcast(xk_tmp(1, 1), 3 * nkstot, mpi_double, 0, inter_pool_comm, ierr)
!
!! convert to fractional coordinates
!do ik = 1, nkstot
!  xk_tmp(:, ik) = matmul(bg_inv, xk_tmp(:, ik))
!end do

!allocate(nk_loc(0:npool-1))
!nk_loc = 0
!nk_loc(rank) = nks
!call mp_sum(nk_loc, inter_pool_comm)
!if (nspin.eq.2) then
!  nk_loc(:) = nk_loc(:)
!endif

!if (nspin.eq.2) then
!  num_kp = nkstot / 2
!else
!  num_kp = nkstot
!endif

!allocate(xk_tmp(3, num_kpoints))
!do i = 1, num_kpoints
!  xk_tmp(:, i) =  matmul(bg_inv, kpoints(:, i))
!enddo
!deallocate(xk_tmp)

! create k-point set
! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice
CALL sirius_create_kset(sctx, num_kpoints, kpoints, wkpoints, .FALSE., ks_handler)

! create ground-state class
CALL sirius_create_ground_state(ks_handler, gs_handler)

END SUBROUTINE setup_sirius


SUBROUTINE update_sirius
USE cell_base, ONLY : alat, at, bg
USE ions_base, ONLY : tau, nat
IMPLICIT NONE
REAL(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp
INTEGER ia
! set lattice vectors of the unit cell (length is in [a.u.])
a1(:) = at(:, 1) * alat
a2(:) = at(:, 2) * alat
a3(:) = at(:, 3) * alat
CALL sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))
!
vlat(:, 1) = a1(:)
vlat(:, 2) = a2(:)
vlat(:, 3) = a3(:)
! get the inverse of Bravais lattice vectors
CALL invert_mtrx(vlat, vlat_inv)
! get the inverse of reciprocal lattice vectors
CALL invert_mtrx(bg, bg_inv)

DO ia = 1, nat
  v1(:) = tau(:, ia) * alat
  ! fractional coordinates
  v1(:) = MATMUL(vlat_inv, v1)
  CALL sirius_set_atom_position(sctx, ia, v1(1))
ENDDO
CALL sirius_update_ground_state(gs_handler)

END SUBROUTINE update_sirius


SUBROUTINE clear_sirius
USE ions_base, ONLY : nsp
IMPLICIT NONE
INTEGER iat

CALL sirius_free_handler(gs_handler)
CALL sirius_free_handler(ks_handler)
CALL sirius_free_handler(sctx)

IF (ALLOCATED(atom_type)) THEN
  DO iat = 1, nsp
    IF (ALLOCATED(atom_type(iat)%qpw)) DEALLOCATE(atom_type(iat)%qpw)
  ENDDO
  DEALLOCATE(atom_type)
ENDIF

END SUBROUTINE


SUBROUTINE get_q_operator_from_sirius
USE uspp_param, ONLY : upf, nh, nhm
USE ions_base,  ONLY : nsp
USE gvect,      ONLY : ngm, mill
USE uspp,       ONLY : qq_nt
IMPLICIT NONE
INTEGER iat, ih, jh, ijh, i

DO iat = 1, nsp
  IF (nh(iat).NE.atom_type(iat)%num_beta_projectors) THEN
     WRITE(*,*) nh(iat), atom_type(iat)%num_beta_projectors
     STOP 'wrong number of beta projectors'
  ENDIF
  IF (upf(iat)%tvanp) THEN
    i = atom_type(iat)%num_beta_projectors
    IF (ALLOCATED(atom_type(iat)%qpw)) DEALLOCATE(atom_type(iat)%qpw)
    ALLOCATE(atom_type(iat)%qpw(ngm, i * (i + 1) / 2))
    ijh = 0
    DO ih = 1, atom_type(iat)%num_beta_projectors
      DO jh = ih, atom_type(iat)%num_beta_projectors
        ijh = ijh + 1
        CALL sirius_get_q_operator(sctx, atom_type(iat)%label, ih, jh, ngm, mill, atom_type(iat)%qpw(:, ijh))
      ENDDO
    ENDDO
  ENDIF
ENDDO

CALL get_q_operator_matrix_from_sirius

END SUBROUTINE get_q_operator_from_sirius


SUBROUTINE get_q_operator_matrix_from_sirius
USE uspp_param, ONLY : upf, nh, nhm
USE ions_base,  ONLY : nsp, ityp, nat
USE uspp,       ONLY : qq_nt, qq_at
IMPLICIT NONE
INTEGER iat, ih, jh, ijh, ia

qq_nt = 0
DO iat = 1, nsp
  IF (nh(iat).NE.atom_type(iat)%num_beta_projectors) THEN
    STOP 'wrong number of beta projectors'
  ENDIF
  IF (upf(iat)%tvanp) THEN
    ijh = 0
    DO ih = 1, atom_type(iat)%num_beta_projectors
      DO jh = ih, atom_type(iat)%num_beta_projectors
        ijh = ijh + 1
        CALL sirius_get_q_operator_matrix(sctx, atom_type(iat)%label, qq_nt(1, 1, iat), nhm)
      ENDDO
    ENDDO
  ENDIF
ENDDO
DO ia = 1, nat
   qq_at(:, :, ia) = qq_nt(:, :, ityp(ia))
END DO

END SUBROUTINE get_q_operator_matrix_from_sirius


SUBROUTINE invert_mtrx(vlat, vlat_inv)
  IMPLICIT NONE
  REAL(8), INTENT(in) :: vlat(3,3)
  REAL(8), INTENT(out) :: vlat_inv(3, 3)
  REAL(8) d1

  d1 = vlat(1,2)*vlat(2,3)*vlat(3,1)-vlat(1,3)*vlat(2,2)*vlat(3,1)+vlat(1,3)*vlat(2,1)*vlat(3,2) &
  &   -vlat(1,1)*vlat(2,3)*vlat(3,2)+vlat(1,1)*vlat(2,2)*vlat(3,3)-vlat(1,2)*vlat(2,1)*vlat(3,3)
  d1 = 1.d0 / d1
  vlat_inv(1,1)=(vlat(2,2)*vlat(3,3)-vlat(2,3)*vlat(3,2))*d1
  vlat_inv(1,2)=(vlat(1,3)*vlat(3,2)-vlat(1,2)*vlat(3,3))*d1
  vlat_inv(1,3)=(vlat(1,2)*vlat(2,3)-vlat(1,3)*vlat(2,2))*d1
  vlat_inv(2,1)=(vlat(2,3)*vlat(3,1)-vlat(2,1)*vlat(3,3))*d1
  vlat_inv(2,2)=(vlat(1,1)*vlat(3,3)-vlat(1,3)*vlat(3,1))*d1
  vlat_inv(2,3)=(vlat(1,3)*vlat(2,1)-vlat(1,1)*vlat(2,3))*d1
  vlat_inv(3,1)=(vlat(2,1)*vlat(3,2)-vlat(2,2)*vlat(3,1))*d1
  vlat_inv(3,2)=(vlat(1,2)*vlat(3,1)-vlat(1,1)*vlat(3,2))*d1
  vlat_inv(3,3)=(vlat(1,1)*vlat(2,2)-vlat(1,2)*vlat(2,1))*d1
END SUBROUTINE


SUBROUTINE get_band_energies_from_sirius
  !
  USE wvfct,    ONLY : nbnd, et
  USE klist,    ONLY : nkstot, nks
  USE lsda_mod, ONLY : nspin
  USE sirius
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !
  REAL(8), ALLOCATABLE :: band_e(:,:)
  INTEGER :: ik, nk, nb, nfv

  ALLOCATE(band_e(nbnd, nkstot))

  ! get band energies
  IF (nspin.NE.2) THEN
    ! non-magnetic or non-collinear case
    DO ik = 1, nkstot
      CALL sirius_get_band_energies(ks_handler, ik, 0, band_e(1, ik))
    END DO
  ELSE
    ! collinear magnetic case
    nk = nkstot / 2
    ! get band energies
    DO ik = 1, nk
      CALL sirius_get_band_energies(ks_handler, ik, 0, band_e(1, ik))
      CALL sirius_get_band_energies(ks_handler, ik, 1, band_e(1, nk + ik))
    END DO

  ENDIF

  ! convert to Ry
  DO ik = 1, nks
    et(:, ik) = 2.d0 * band_e(:, global_kpoint_index(nkstot, ik))
  ENDDO

  DEALLOCATE(band_e)

END SUBROUTINE get_band_energies_from_sirius


SUBROUTINE put_band_occupancies_to_sirius
  !
  USE wvfct,    ONLY : nbnd, wg
  USE klist,    ONLY : nkstot, nks, wk
  USE lsda_mod, ONLY : nspin
  USE mp_pools, ONLY : inter_pool_comm
  USE parallel_include
  USE sirius
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !
  REAL(8), ALLOCATABLE :: bnd_occ(:, :)
  REAL(8) :: maxocc
  INTEGER :: ik, ierr, nk, nb

  ! compute occupancies
  ALLOCATE(bnd_occ(nbnd, nkstot))
  bnd_occ = 0.d0
  ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
  maxocc = 2.d0
  IF (nspin.GT.1) THEN
    maxocc = 1.d0
  ENDIF
  DO ik = 1, nks
    bnd_occ(:, global_kpoint_index(nkstot, ik)) = maxocc * wg(:, ik) / wk(ik)
  ENDDO
  CALL mpi_allreduce(MPI_IN_PLACE, bnd_occ(1, 1), nbnd * nkstot, MPI_DOUBLE, MPI_SUM, inter_pool_comm, ierr)

  IF (nspin.NE.2) THEN
    ! set band occupancies
    DO ik = 1, nkstot
      CALL sirius_set_band_occupancies(ks_handler, ik, 0, bnd_occ(1, ik))
    ENDDO
  ELSE
    nk = nkstot / 2
    DO ik = 1, nk
      CALL sirius_set_band_occupancies(ks_handler, ik, 0, bnd_occ(1, ik))
      CALL sirius_set_band_occupancies(ks_handler, ik, 1, bnd_occ(1, ik + nk))
    ENDDO
  ENDIF

  DEALLOCATE(bnd_occ)

END SUBROUTINE put_band_occupancies_to_sirius


SUBROUTINE get_band_occupancies_from_sirius
  !
  USE wvfct,    ONLY : nbnd, wg
  USE klist,    ONLY : nkstot, nks, wk
  USE lsda_mod, ONLY : nspin
  USE mp_pools, ONLY : inter_pool_comm
  USE parallel_include
  USE sirius
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: global_kpoint_index
  !
  REAL(8), ALLOCATABLE :: bnd_occ(:, :)
  REAL(8) :: maxocc
  INTEGER :: ik, ierr, nk, nb

  ! compute occupancies
  ALLOCATE(bnd_occ(nbnd, nkstot))
  IF (nspin.NE.2) THEN
    ! set band occupancies
    DO ik = 1, nkstot
      CALL sirius_get_band_occupancies(ks_handler, ik, 0, bnd_occ(1, ik))
    ENDDO
  ELSE
    nk = nkstot / 2
    DO ik = 1, nk
      CALL sirius_get_band_occupancies(ks_handler, ik, 0, bnd_occ(1, ik))
      CALL sirius_get_band_occupancies(ks_handler, ik, 1, bnd_occ(1, ik + nk))
    ENDDO
  ENDIF

  ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
  maxocc = 2.d0
  IF (nspin.GT.1) THEN
    maxocc = 1.d0
  ENDIF
  DO ik = 1, nks
    wg(:, ik) = bnd_occ(:, global_kpoint_index(nkstot, ik)) / maxocc * wk(ik)
  ENDDO

  DEALLOCATE(bnd_occ)

END SUBROUTINE get_band_occupancies_from_sirius



SUBROUTINE get_density_matrix_from_sirius
  !
  USE scf,        ONLY : rho
  USE ions_base,  ONLY : nat, nsp, ityp
  USE uspp_param, ONLY : nhm, nh
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE
  !
  COMPLEX(8), ALLOCATABLE :: dens_mtrx(:,:,:)
  INTEGER iat, na, ijh, ih, jh, ispn
  ! complex density matrix in SIRIUS has at maximum three components
  ALLOCATE(dens_mtrx(nhm, nhm, 3))
  DO iat = 1, nsp
    DO na = 1, nat
      IF (ityp(na).EQ.iat.AND.ALLOCATED(rho%bec)) THEN
        rho%bec(:, na, :) = 0.d0
        CALL sirius_get_density_matrix(gs_handler, na, dens_mtrx(1, 1, 1), nhm)
        ijh = 0
        DO ih = 1, nh(iat)
          DO jh = ih, nh(iat)
            ijh = ijh + 1
            IF (nspin.LE.2) THEN
              DO ispn = 1, nspin
                rho%bec(ijh, na, ispn) = dreal(dens_mtrx(ih, jh, ispn))
              ENDDO
            ENDIF
            IF (nspin.EQ.4) THEN
              rho%bec(ijh, na, 1) = dreal(dens_mtrx(ih, jh, 1) + dens_mtrx(ih, jh, 2))
              rho%bec(ijh, na, 4) = dreal(dens_mtrx(ih, jh, 1) - dens_mtrx(ih, jh, 2))
              rho%bec(ijh, na, 2) = 2.d0 * dreal(dens_mtrx(ih, jh, 3))
              rho%bec(ijh, na, 3) = -2.d0 * dimag(dens_mtrx(ih, jh, 3))
            ENDIF
            ! off-diagonal elements have a weight of 2
            IF (ih.NE.jh) THEN
              DO ispn = 1, nspin
                rho%bec(ijh, na, ispn) = rho%bec(ijh, na, ispn) * 2.d0
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(dens_mtrx)
END SUBROUTINE get_density_matrix_from_sirius


SUBROUTINE get_density_from_sirius
  !
  USE scf,        ONLY : rho
  USE gvect,      ONLY : mill, ngm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE lsda_mod,   ONLY : nspin
  USE ions_base,  ONLY : nat, nsp, ityp
  USE uspp_param, ONLY : nhm, nh
  USE sirius
  !
  IMPLICIT NONE
  !
  INTEGER iat, ig, ih, jh, ijh, na, ispn
  COMPLEX(8) z1, z2

  ! get rho(G)
  CALL sirius_get_pw_coeffs(gs_handler, "rho", rho%of_g(:, 1), ngm, mill, intra_bgrp_comm)
  IF (nspin.EQ.2) THEN
    CALL sirius_get_pw_coeffs(gs_handler, "magz", rho%of_g(:, 2), ngm, mill, intra_bgrp_comm)
    !! convert to rho_{up}, rho_{dn}
    !do ig = 1, ngm
    !  z1 = rho%of_g(ig, 1)
    !  z2 = rho%of_g(ig, 2)
    !  rho%of_g(ig, 1) = 0.5 * (z1 + z2)
    !  rho%of_g(ig, 2) = 0.5 * (z1 - z2)
    !enddo
  ENDIF
  IF (nspin.EQ.4) THEN
    CALL sirius_get_pw_coeffs(gs_handler, "magx", rho%of_g(:, 2), ngm, mill, intra_bgrp_comm)
    CALL sirius_get_pw_coeffs(gs_handler, "magy", rho%of_g(:, 3), ngm, mill, intra_bgrp_comm)
    CALL sirius_get_pw_coeffs(gs_handler, "magz", rho%of_g(:, 4), ngm, mill, intra_bgrp_comm)
  ENDIF
  ! get density matrix
  CALL get_density_matrix_from_sirius
END SUBROUTINE get_density_from_sirius


SUBROUTINE put_density_to_sirius
  !
  USE scf,        ONLY : rho
  USE gvect,      ONLY : mill, ngm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE lsda_mod,   ONLY : nspin
  USE ions_base,  ONLY : nat, nsp, ityp
  USE uspp_param, ONLY : nhm, nh
  USE sirius
  IMPLICIT NONE
  !
  COMPLEX(8), ALLOCATABLE :: rho_tot(:), mag(:)
  INTEGER iat, ig, ih, jh, ijh, na, ispn
  REAL(8) :: fact
  !
  !if (nspin.eq.1.or.nspin.eq.4) then
    CALL sirius_set_pw_coeffs(gs_handler, "rho", rho%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm)
  !endif

  IF (nspin.EQ.2) THEN
    !allocate(rho_tot(ngm))
    !allocate(mag(ngm))
    !do ig = 1, ngm
    !  rho_tot(ig) = rho%of_g(ig, 1) + rho%of_g(ig, 2)
    !  mag(ig) = rho%of_g(ig, 1) - rho%of_g(ig, 2)
    !enddo
    !call sirius_set_pw_coeffs(gs_handler, "rho", rho_tot(1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    !call sirius_set_pw_coeffs(gs_handler, "magz", mag(1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    !deallocate(rho_tot)
    !deallocate(mag)
    CALL sirius_set_pw_coeffs(gs_handler, "magz", rho%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF

  IF (nspin.EQ.4) THEN
    CALL sirius_set_pw_coeffs(gs_handler, "magx", rho%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm)
    CALL sirius_set_pw_coeffs(gs_handler, "magy", rho%of_g(:, 3), .TRUE., ngm, mill, intra_bgrp_comm)
    CALL sirius_set_pw_coeffs(gs_handler, "magz", rho%of_g(:, 4), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF
END SUBROUTINE put_density_to_sirius


SUBROUTINE put_density_matrix_to_sirius
  !
  USE scf,        ONLY : rho
  USE ions_base,  ONLY : nat, nsp, ityp
  USE lsda_mod,   ONLY : nspin
  USE uspp_param, ONLY : nhm, nh
  USE uspp,       ONLY : becsum
  IMPLICIT NONE
  !
  INTEGER iat, na, ih, jh, ijh, ispn
  COMPLEX(8), ALLOCATABLE :: dens_mtrx(:,:,:)
  REAL(8), ALLOCATABLE :: dens_mtrx_tmp(:, :, :)
  REAL(8) fact
  ! set density matrix
  ! complex density matrix in SIRIUS has at maximum three components
  ALLOCATE(dens_mtrx_tmp(nhm * (nhm + 1) / 2, nat, nspin))
  !if (allocated(rho%bec)) then
  !  dens_mtrx_tmp = rho%bec
  !else
    dens_mtrx_tmp = becsum
  !endif

  ALLOCATE(dens_mtrx(nhm, nhm, 3))
  DO iat = 1, nsp
    DO na = 1, nat
      IF (ityp(na).EQ.iat) THEN
        dens_mtrx = (0.d0, 0.d0)
        ijh = 0
        DO ih = 1, nh(iat)
          DO jh = ih, nh(iat)
            ijh = ijh + 1
            ! off-diagonal elements have a weight of 2
            IF (ih.NE.jh) THEN
              fact = 0.5d0
            ELSE
              fact = 1.d0
            ENDIF
            IF (nspin.LE.2) THEN
              DO ispn = 1, nspin
                dens_mtrx(ih, jh, ispn) = fact * dens_mtrx_tmp(ijh, na, ispn)
                dens_mtrx(jh, ih, ispn) = fact * dens_mtrx_tmp(ijh, na, ispn)
              ENDDO
            ENDIF
            IF (nspin.EQ.4) THEN
              ! 0.5 * (rho + mz)
              dens_mtrx(ih, jh, 1) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) + dens_mtrx_tmp(ijh, na, 4))
              dens_mtrx(jh, ih, 1) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) + dens_mtrx_tmp(ijh, na, 4))
              ! 0.5 * (rho - mz)
              dens_mtrx(ih, jh, 2) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) - dens_mtrx_tmp(ijh, na, 4))
              dens_mtrx(jh, ih, 2) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) - dens_mtrx_tmp(ijh, na, 4))
              ! 0.5 * (mx - I * my)
              dens_mtrx(ih, jh, 3) = fact * 0.5 * dcmplx(dens_mtrx_tmp(ijh, na, 2), -dens_mtrx_tmp(ijh, na, 3))
              dens_mtrx(jh, ih, 3) = fact * 0.5 * dcmplx(dens_mtrx_tmp(ijh, na, 2), -dens_mtrx_tmp(ijh, na, 3))
            ENDIF
          ENDDO
        ENDDO
        CALL sirius_set_density_matrix(gs_handler, na, dens_mtrx(1, 1, 1), nhm)
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(dens_mtrx)
  DEALLOCATE(dens_mtrx_tmp)
END SUBROUTINE put_density_matrix_to_sirius


SUBROUTINE put_d_matrix_to_sirius
USE uspp_param,           ONLY : nhm
USE ions_base,            ONLY : nat
USE lsda_mod,             ONLY : nspin
USE uspp,                 ONLY : deeq
IMPLICIT NONE
REAL(8), ALLOCATABLE :: deeq_tmp(:,:)
INTEGER ia, is
ALLOCATE(deeq_tmp(nhm, nhm))
DO ia = 1, nat
  DO is = 1, nspin
    IF (nspin.EQ.2.AND.is.EQ.1) THEN
      deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) + deeq(:, :, ia, 2)) / 2 ! convert to Ha
    ENDIF
    IF (nspin.EQ.2.AND.is.EQ.2) THEN
      deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) - deeq(:, :, ia, 2)) / 2 ! convert to Ha
    ENDIF
    IF (nspin.EQ.1.OR.nspin.EQ.4) THEN
      deeq_tmp(:, :) = deeq(:, :, ia, is) / 2 ! convert to Ha
    ENDIF
    CALL sirius_set_d_operator_matrix(sctx, ia, is, deeq_tmp(1, 1), nhm)
  ENDDO
ENDDO
DEALLOCATE(deeq_tmp)
END SUBROUTINE put_d_matrix_to_sirius


SUBROUTINE get_d_matrix_from_sirius
USE uspp_param,           ONLY : nhm
USE ions_base,            ONLY : nat
USE lsda_mod,             ONLY : nspin
USE uspp,                 ONLY : deeq
IMPLICIT NONE
REAL(8) d1, d2
INTEGER ia, is, i, j
! get D-operator matrix
DO ia = 1, nat
  DO is = 1, nspin
    CALL sirius_get_d_operator_matrix(sctx, ia, is, deeq(1, 1, ia, is), nhm)
  ENDDO
  IF (nspin.EQ.2) THEN
    DO i = 1, nhm
      DO j = 1, nhm
        d1 = deeq(i, j, ia, 1)
        d2 = deeq(i, j, ia, 2)
        deeq(i, j, ia, 1) = d1 + d2
        deeq(i, j, ia, 2) = d1 - d2
      ENDDO
    ENDDO
  ENDIF
  ! convert to Ry
  deeq(:, :, ia, :) = deeq(:, :, ia, :) * 2
ENDDO
END SUBROUTINE get_d_matrix_from_sirius

SUBROUTINE put_potential_to_sirius
  USE scf,                  ONLY : v, vltot
  USE gvect,                ONLY : mill, ngm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : nspin_mag
  USE wavefunctions,        ONLY : psic

  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE uspp,                 ONLY : deeq
  USE uspp_param,           ONLY : nhm
  USE paw_variables,        ONLY : okpaw
  USE ions_base,            ONLY : nat
  USE sirius
  !
  IMPLICIT NONE
  !
  COMPLEX(8), ALLOCATABLE :: vxcg(:)
  COMPLEX(8) :: z1, z2
  INTEGER ig, is, ir, i, ia, j
  CHARACTER(10) label
  REAL(8), ALLOCATABLE :: deeq_tmp(:,:)
  REAL(8) :: d1,d2
  !
  IF (nspin.EQ.1.OR.nspin.EQ.4) THEN
    ! add local part of the potential and transform to PW domain
    psic(:) = v%of_r(:, 1) + vltot(:)
    CALL fwfft('Rho', psic, dfftp)
    ! convert to Hartree
    DO ig = 1, ngm
      v%of_g(ig, 1) = psic(dfftp%nl(ig)) * 0.5d0
    ENDDO
    ! set effective potential
    CALL sirius_set_pw_coeffs(gs_handler, "veff", v%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF

  IF (nspin.EQ.2) THEN
    DO is = 1, 2
      ! add local part of the potential and transform to PW domain
      psic(:) = v%of_r(:, is) + vltot(:)
      CALL fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      DO ig = 1, ngm
         v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      ENDDO
    ENDDO

    DO ig = 1, ngm
      z1 = v%of_g(ig, 1)
      z2 = v%of_g(ig, 2)
      v%of_g(ig, 1) = 0.5 * (z1 + z2)
      v%of_g(ig, 2) = 0.5 * (z1 - z2)
    ENDDO
    ! set effective potential and magnetization
    CALL sirius_set_pw_coeffs(gs_handler, "veff", v%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm)
    CALL sirius_set_pw_coeffs(gs_handler, "bz",   v%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF

  IF (nspin.EQ.4) THEN
    DO is = 2, nspin_mag
      psic(:) = v%of_r(:, is)
      CALL fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      DO ig = 1, ngm
        v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      ENDDO
      IF (is.EQ.2) label="bx"
      IF (is.EQ.3) label="by"
      IF (is.EQ.4) label="bz"

      CALL sirius_set_pw_coeffs(gs_handler, label, v%of_g(:, is), .TRUE., ngm, mill, intra_bgrp_comm)
    ENDDO
  ENDIF

!  ! convert Vxc to plane-wave domain
!  if (nspin.eq.1.or.nspin.eq.4) then
!    do ir = 1, dfftp%nnr
!      psic(ir) = vxc(ir, 1)
!    enddo
!  else
!    do ir = 1, dfftp%nnr
!      psic(ir) = 0.5d0 * (vxc(ir, 1) + vxc(ir, 2))
!    enddo
!  endif
!  call fwfft('Rho', psic, dfftp)
!  allocate(vxcg(ngm))
!  ! convert to Hartree
!  do ig = 1, ngm
!     vxcg(ig) = psic(dfftp%nl(ig)) * 0.5d0
!  end do
!  ! set XC potential
!  call sirius_set_pw_coeffs("vxc", vxcg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vxcg)
!
!  ! update D-operator matrix
!  !call sirius_generate_d_operator_matrix()
!  !if (okpaw) then
!    allocate(deeq_tmp(nhm, nhm))
!  !  !! get D-operator matrix
!  !  !do ia = 1, nat
!  !  !  do is = 1, nspin
!  !  !    call sirius_get_d_operator_matrix(ia, is, deeq(1, 1, ia, is), nhm)
!  !  !  enddo
!  !  !  if (nspin.eq.2) then
!  !  !    do i = 1, nhm
!  !  !      do j = 1, nhm
!  !  !        d1 = deeq(i, j, ia, 1)
!  !  !        d2 = deeq(i, j, ia, 2)
!  !  !        deeq(i, j, ia, 1) = d1 + d2
!  !  !        deeq(i, j, ia, 2) = d1 - d2
!  !  !      enddo
!  !  !    enddo
!  !  !  endif
!  !  !  ! convert to Ry
!  !  !  deeq(:, :, ia, :) = deeq(:, :, ia, :) * 2
!  !  !enddo
!  !  !call add_paw_to_deeq(deeq)
!    do ia = 1, nat
!      do is = 1, nspin
!        if (nspin.eq.2.and.is.eq.1) then
!          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) + deeq(:, :, ia, 2)) / 2 ! convert to Ha
!        endif
!        if (nspin.eq.2.and.is.eq.2) then
!          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) - deeq(:, :, ia, 2)) / 2 ! convert to Ha
!        endif
!        if (nspin.eq.1.or.nspin.eq.4) then
!          deeq_tmp(:, :) = deeq(:, :, ia, is) / 2 ! convert to Ha
!        endif
!        call sirius_set_d_operator_matrix(ia, is, deeq_tmp(1, 1), nhm)
!      enddo
!    enddo
!    deallocate(deeq_tmp)
!  !endif

END SUBROUTINE put_potential_to_sirius

SUBROUTINE put_q_operator_matrix_to_sirius
USE uspp,       ONLY : qq_nt
USE uspp_param, ONLY : upf, nhm
USE ions_base,  ONLY : nsp, atm
IMPLICIT NONE
INTEGER iat

DO iat = 1, nsp
  IF (upf(iat)%tvanp) THEN
    CALL sirius_set_q_operator_matrix(sctx, atom_type(iat)%label, qq_nt(1, 1, iat), nhm)
  ENDIF
ENDDO

END SUBROUTINE put_q_operator_matrix_to_sirius


SUBROUTINE get_wave_functions_from_sirius
USE klist, ONLY : nkstot, nks, ngk, igk_k
USE gvect, ONLY : mill
USE buffers, ONLY : save_buffer
USE io_files, ONLY : iunwfc, nwordwfc
USE bp, ONLY : lelfield
USE noncollin_module, ONLY : npol
USE wvfct, ONLY : npwx, nbnd
USE wavefunctions, ONLY : evc
USE lsda_mod, ONLY : isk, lsda
USE mp_pools, ONLY : inter_pool_comm
USE parallel_include
IMPLICIT NONE
INTEGER, EXTERNAL :: global_kpoint_index
INTEGER, ALLOCATABLE :: gvl(:,:)
INTEGER ig, ik, ik_, i, j, ispn, rank, ierr, nksmax
COMPLEX(8) z1

! rank of communicator that distributes k-points
CALL mpi_comm_rank(inter_pool_comm, rank, ierr)
CALL mpi_allreduce(nks, nksmax, 1, MPI_INTEGER, MPI_MAX, inter_pool_comm, ierr)

!CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )

ALLOCATE(gvl(3, npwx))
DO ik = 1, nksmax
  IF (ik.LE.nks) THEN
    DO ig = 1, ngk(ik)
      gvl(:,ig) = mill(:, igk_k(ig, ik))
    ENDDO
    !
    ik_ = global_kpoint_index(nkstot, ik)
    ispn = isk(ik)
    IF (lsda.AND.ispn.EQ.2) THEN
      ik_ = ik_ - nkstot / 2
    ENDIF
    CALL sirius_get_wave_functions(ks_handler, ik_, ispn, ngk(ik), gvl(1, 1), evc(1, 1), npwx, npol)
    IF (nks > 1 .OR. lelfield) THEN
      CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
    ENDIF
  ELSE
    CALL sirius_get_wave_functions(ks_handler, -1, -1, -1, -1, z1, -1, -1)
  ENDIF
ENDDO
DEALLOCATE(gvl)

END SUBROUTINE get_wave_functions_from_sirius


!subroutine put_vltot_to_sirius
!  use scf,       only : vltot
!  use gvect, only : mill, ngm
!  use wavefunctions_module, only : psic
!  use fft_interfaces,       only : fwfft, invfft
!  use fft_base,             only : dfftp
!  use mp_bands, only : intra_bgrp_comm
!  use sirius
!  !
!  implicit none
!  !
!  complex(8), allocatable :: vg(:)
!  integer ig
!  !
!  allocate(vg(ngm))
!  psic(:) = vltot(:)
!  call fwfft('Rho', psic, dfftp)
!  ! convert to Hartree
!  do ig = 1, ngm
!    vg(ig) = psic(dfftp%nl(ig)) * 0.5d0 ! convert to Ha
!  enddo
!  ! set local potential
!  call sirius_set_pw_coeffs("vloc", vg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vg)
!end subroutine put_vltot_to_sirius

!subroutine get_rhoc_from_sirius
!  use uspp_param,only : upf
!  use ener,      only : etxcc
!  use scf,       only : rho_core, rhog_core
!  use control_flags, only : gamma_only
!  use wavefunctions_module, only : psic
!  use gvect, only : mill, ngm
!  use scf, only : rho_core, rhog_core
!  use mp_bands, only : intra_bgrp_comm
!  use ions_base, only : ntyp => nsp
!  use fft_interfaces,only : invfft
!  use fft_base,  only : dfftp
!  use sirius
!  !
!  implicit none
!
!  etxcc = 0.0d0
!  if (any(upf(1:ntyp)%nlcc)) then
!    call sirius_get_pw_coeffs("rhoc", rhog_core(1), ngm, mill(1, 1), intra_bgrp_comm)
!    psic(:) = (0.d0, 0.d0)
!    psic(dfftp%nl(:)) = rhog_core(:)
!    if (gamma_only) psic(dfftp%nlm(:)) = conjg(rhog_core(:))
!    call invfft ('Rho', psic, dfftp)
!    rho_core(:) = psic(:)
!  else
!    rhog_core(:) = 0.0d0
!    rho_core(:)  = 0.0d0
!  endif
!
!end subroutine get_rhoc_from_sirius

!subroutine set_vloc_sirius
!use sirius
!use gvect, only : ngm, mill, igtongl, ngl
!use vlocal, only : vloc
!use mp_bands, only : intra_bgrp_comm
!use ions_base, only : atm
!integer nt,i
!real(8), allocatable :: tmp(:)
!
!allocate(tmp(ngm))
!vloc(:,:) = 0.d0
!do nt = 1, ntyp
!  call sirius_get_pw_coeffs_real(string(atm(nt)), "vloc", tmp(1), ngm, mill(1, 1), intra_bgrp_comm)
!  do i = 1, ngm
!    vloc(igtongl(i), nt) = tmp(i) * 2 ! convert to Ry
!  enddo
!enddo
!
!deallocate(tmp, tmp1)
!
!call set_vloc_sirius
!CALL setlocal()
!end subroutine

!subroutine get_vloc_from_sirius
!  use wavefunctions_module, only : psic
!  use gvect, only : mill, ngm, gg
!  use scf, only: vltot, v_of_0
!  use fft_interfaces, only : fwfft, invfft
!  use fft_base, only : dfftp
!  use constants, only : eps8
!  use control_flags, only : gamma_only
!  use mp_bands, only : intra_bgrp_comm
!  use mp, only : mp_bcast, mp_sum
!  use sirius
!  !
!  implicit none
!  !
!  complex(8), allocatable :: vpw(:)
!  allocate(vpw(ngm))
!  call sirius_get_pw_coeffs("vloc", vpw(1), ngm, mill(1, 1), intra_bgrp_comm)
!  psic(:) = 0.d0
!  psic(dfftp%nl(:)) = vpw(:)
!  if (gamma_only) psic(dfftp%nlm(:)) = conjg(vpw(:))
!  call invfft('Rho', psic, dfftp)
!  vltot(:) = dble(psic(:)) * 2 ! convert to Ry
!  v_of_0=0.d0
!  IF (gg(1) < eps8) v_of_0 = dble(vpw(1))
!  !
!  call mp_sum(v_of_0, intra_bgrp_comm)
!  deallocate(vpw)
!
!end subroutine get_vloc_from_sirius

!subroutine get_density_matrix_from_sirius
!implicit none
!real(8), allocatable :: dens_mtrx(:, :, :)
!integer iat, na
!
!allocate(dens_mtrx(nhm, nhm, 3))
!do iat = 1, nsp
!  do na = 1, nat
!    if (ityp(na).eq.iat.and.allocated(rho%bec)) then
!      rho%bec(:, na, :) = 0.d0
!      call sirius_get_density_matrix(na, dens_mtrx(1, 1, 1), nhm)
!
!      ijh = 0
!      do ih = 1, nh(iat)
!        do jh = ih, nh(iat)
!          ijh = ijh + 1
!          if (nspin.le.2) then
!            do ispn = 1, nspin
!              rho%bec(ijh, na, ispn) = dreal(dens_mtrx(ih, jh, ispn))
!            enddo
!          endif
!          if (nspin.eq.4) then
!            rho%bec(ijh, na, 1) = dreal(dens_mtrx(ih, jh, 1) + dens_mtrx(ih, jh, 2))
!            rho%bec(ijh, na, 4) = dreal(dens_mtrx(ih, jh, 1) - dens_mtrx(ih, jh, 2))
!            rho%bec(ijh, na, 2) = 2.d0 * dreal(dens_mtrx(ih, jh, 3))
!            rho%bec(ijh, na, 3) = -2.d0 * dimag(dens_mtrx(ih, jh, 3))
!          endif
!          ! off-diagonal elements have a weight of 2
!          if (ih.ne.jh) then
!            do ispn = 1, nspin
!              rho%bec(ijh, na, ispn) = rho%bec(ijh, na, ispn) * 2.d0
!            enddo
!          endif
!        enddo
!      enddo
!    endif
!  enddo
!enddo
!deallocate(dens_mtrx)
!
!end subroutine get_density_matrix_from_sirius

SUBROUTINE put_xc_functional_to_sirius
IMPLICIT NONE

  IF (get_meta().NE.0.OR.get_inlc().NE.0) THEN
   WRITE(*,*)get_igcx()
   WRITE(*,*)get_igcc()
   WRITE(*,*)get_meta()
   WRITE(*,*)get_inlc()
   STOP ("interface for this XC functional is not implemented")
  ENDIF

  IF (get_iexch().NE.0.AND.get_igcx().EQ.0) THEN
   WRITE(*,*) 'iexch, igcx:', get_iexch(), get_igcx()
   SELECT CASE(get_iexch())
   CASE(0)
   CASE(1)
     CALL sirius_add_xc_functional(sctx, "XC_LDA_X")
   CASE default
     STOP ("interface for this exchange functional is not implemented")
   END SELECT
  ENDIF

  IF (get_iexch().NE.0.AND.get_igcx().NE.0) THEN
   SELECT CASE(get_igcx())
   CASE(0)
   CASE(2)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PW91")
   CASE(3)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE")
   CASE(10)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE_SOL")
   CASE default
     WRITE(*,*)get_igcx()
     STOP ("interface for this gradient exchange functional is not implemented")
   END SELECT
  ENDIF

  IF (get_icorr().NE.0.AND.get_igcc().EQ.0) THEN
   SELECT CASE(get_icorr())
   CASE(0)
   CASE(1)
     CALL sirius_add_xc_functional(sctx, "XC_LDA_C_PZ")
   CASE(4)
     CALL sirius_add_xc_functional(sctx, "XC_LDA_C_PW")
   CASE default
     STOP ("interface for this correlation functional is not implemented")
   END SELECT
  ENDIF

  IF (get_icorr().NE.0.AND.get_igcc().NE.0) THEN
   SELECT CASE(get_igcc())
   CASE(0)
   CASE(2)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PW91")
   CASE(4)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PBE")
   CASE(8)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PBE_SOL")
   CASE default
     STOP ("interface for this gradient correlation functional is not implemented")
   END SELECT
  ENDIF
END SUBROUTINE put_xc_functional_to_sirius


SUBROUTINE insert_xc_functional_to_sirius
  IMPLICIT NONE

  IF (get_meta().NE.0.OR.get_inlc().NE.0) THEN
    WRITE(*,*)get_igcx()
    WRITE(*,*)get_igcc()
    WRITE(*,*)get_meta()
    WRITE(*,*)get_inlc()
    STOP ("interface for this XC functional is not implemented")
  ENDIF

  IF (get_iexch().NE.0.AND.get_igcx().EQ.0) THEN
    SELECT CASE(get_iexch())
    CASE(0)
    CASE(1)
      CALL sirius_insert_xc_functional(gs_handler, "XC_LDA_X")
    CASE default
      STOP ("interface for this exchange functional is not implemented")
    END SELECT
  ENDIF

  IF (get_iexch().NE.0.AND.get_igcx().NE.0) THEN
    SELECT CASE(get_igcx())
    CASE(0)
    CASE(2)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_X_PW91")
    CASE(3)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_X_PBE")
    CASE(10)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_X_PBE_SOL")
    CASE(21)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_X_PW86")
    CASE(22)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_X_B86")
    CASE default
      WRITE(*,*)get_igcx()
      STOP ("interface for this gradient exchange functional is not implemented")
    END SELECT
  ENDIF

  IF (get_icorr().NE.0.AND.get_igcc().EQ.0) THEN
    SELECT CASE(get_icorr())
    CASE(0)
    CASE(1)
      CALL sirius_insert_xc_functional(gs_handler, "XC_LDA_C_PZ")
    CASE(4)
      CALL sirius_insert_xc_functional(gs_handler, "XC_LDA_C_PW")
    CASE default
      STOP ("interface for this correlation functional is not implemented")
    END SELECT
  ENDIF

  IF (get_icorr().NE.0.AND.get_igcc().NE.0) THEN
    SELECT CASE(get_igcc())
    CASE(0)
    CASE(1)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_C_PW86")
    CASE(2)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_C_PW91")
    CASE(4)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_C_PBE")
    CASE(8)
      CALL sirius_insert_xc_functional(gs_handler, "XC_GGA_C_PBE_SOL")
    CASE default
      STOP ("interface for this gradient correlation functional is not implemented")
    END SELECT
  ENDIF
END SUBROUTINE insert_xc_functional_to_sirius


SUBROUTINE write_json()
USE ener
USE force_mod
USE ions_base
IMPLICIT NONE
INTEGER i

OPEN(200, file="output.json", action="write", form="formatted")
WRITE(200,'("{")')
WRITE(200,'("  ""energy"": {")')
WRITE(200,'("    ""total"": ", G18.10)')etot
WRITE(200, '("  },")')
WRITE(200,'("  ""stress"": {")')
WRITE(200,'("    ""total"": [")')
WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') sigma(1, 1), sigma(1, 2), sigma(1, 3)
WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') sigma(2, 1), sigma(2, 2), sigma(2, 3)
WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "]")')  sigma(3, 1), sigma(3, 2), sigma(3, 3)
WRITE(200,'("    ]")')
WRITE(200,'("  },")')
WRITE(200,'("  ""force"": {")')
WRITE(200,'("    ""total"": [")')
DO i = 1, nat
  IF (i.EQ.nat) THEN
    WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "]")') force(1, i), force(2, i), force(3, i)
  ELSE
    WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') force(1, i), force(2, i), force(3, i)
  ENDIF
ENDDO
WRITE(200,'("    ]")')
WRITE(200,'("  }")')
WRITE(200,'("}")')
CLOSE(200)

END SUBROUTINE

FUNCTION idx_m_qe(m) RESULT(m1)
  IMPLICIT NONE
  INTEGER :: m
  INTEGER :: m1

  IF (m .GT. 0) THEN
     m1 = 2 * m - 1
  ELSE
     m1 = -2 * m
  ENDIF
END FUNCTION idx_m_qe

SUBROUTINE qe_to_sirius_real(ns, ns_sirius)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin

  COMPLEX(8), INTENT(out) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  REAL(8), INTENT(in) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  INTEGER :: m1, m2, mm2, mm1, is, nt, na
  ns_sirius(:, :, :, :) = 0.0
  DO na = 1, nat
     nt = ityp (na)
     IF (is_hubbard(nt)) THEN
        DO m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           DO m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              DO is = 1, nspin
                 ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, is, na) = ns(mm1 + 1, mm2 + 1, is, na)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
END SUBROUTINE qe_to_sirius_real

SUBROUTINE qe_to_sirius_complex(ns, ns_sirius)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin

  COMPLEX(8), INTENT(out) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  COMPLEX(8), INTENT(in) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  INTEGER :: m1, m2, mm2, mm1, is, nt, na

  ns_sirius(:, :, :, :) = 0.0
  DO na = 1, nat
     nt = ityp (na)
     IF (is_hubbard(nt)) THEN
        DO m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           DO m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 1, na) = ns(mm1 + 1, mm2 + 1, 1, na)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 2, na) = ns(mm1 + 1, mm2 + 1, 4, na)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 3, na) = ns(mm1 + 1, mm2 + 1, 2, na)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 4, na) = ns(mm1 + 1, mm2 + 1, 3, na)
           ENDDO
        ENDDO
     ENDIF
  ENDDO
END SUBROUTINE qe_to_sirius_complex

SUBROUTINE sirius_to_qe_real(ns_sirius, ns)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin
  INTEGER :: m1, m2, mm2, mm1, is, nt, na

  COMPLEX(8), INTENT(in) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  REAL(8), INTENT(out) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)

  ns(:, :, :, :) = 0.0
  DO na = 1, nat
     nt = ityp (na)
     IF (is_hubbard(nt)) THEN
        DO m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           DO m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              DO is = 1, nspin
                 ns(mm1 + 1, mm2 + 1, is, na) = REAL(ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, is, na))
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
END SUBROUTINE sirius_to_qe_real

SUBROUTINE sirius_to_qe_complex(ns_sirius, ns)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin

  COMPLEX(8), INTENT(in) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat)
  COMPLEX(8), INTENT(out) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  INTEGER :: m1, m2, mm2, mm1, is, nt, na

  ns(:, :, :, :) = 0.0
  DO na = 1, nat
     nt = ityp (na)
     IF (is_hubbard(nt)) THEN
        DO m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           DO m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              ns(mm1 + 1, mm2 + 1, 1, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 1, na)
              ns(mm1 + 1, mm2 + 1, 4, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 2, na)
              ns(mm1 + 1, mm2 + 1, 2, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 3, na)
              ns(mm1 + 1, mm2 + 1, 3, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 4, na)
           ENDDO
        ENDDO
     ENDIF
  ENDDO
END SUBROUTINE sirius_to_qe_complex

SUBROUTINE qe_sirius_set_hubbard_occupancy(rho)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE lsda_mod,             ONLY : nspin

  IMPLICIT NONE

  TYPE(scf_type), INTENT(IN) :: rho  ! the valence charge
  COMPLEX(8), ALLOCATABLE :: ns_sirius(:, :, :, :)


  IF (noncolin) THEN
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     CALL qe_to_sirius_complex(rho%ns_nc(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  ELSE
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat))
     CALL qe_to_sirius_real(rho%ns(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  ENDIF
  CALL sirius_set_hubbard_occupancies(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)
  DEALLOCATE(ns_sirius)
END SUBROUTINE qe_sirius_set_hubbard_occupancy

SUBROUTINE qe_sirius_get_hubbard_occupancy(rho)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE lsda_mod,             ONLY : nspin
  IMPLICIT NONE

  TYPE(scf_type), INTENT(out) :: rho  ! the valence charge
  COMPLEX(8), ALLOCATABLE :: ns_sirius(:, :, :, :)


  IF (noncolin) THEN
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     CALL sirius_get_hubbard_occupancies(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)
     CALL sirius_to_qe_complex(ns_sirius(1, 1, 1, 1), rho%ns_nc(1, 1, 1, 1))
  ELSE
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 2, nat))
     CALL sirius_get_hubbard_occupancies(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)
     CALL sirius_to_qe_real(ns_sirius(1, 1, 1, 1), rho%ns(1, 1, 1, 1))
  ENDIF

  DEALLOCATE(ns_sirius)
END SUBROUTINE qe_sirius_get_hubbard_occupancy

SUBROUTINE qe_sirius_set_hubbard_potential(v)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  IMPLICIT NONE
  TYPE(scf_type), INTENT(IN) :: v  ! the valence charge
  COMPLEX(8), ALLOCATABLE :: ns_sirius(:, :, :, :)

  IF (noncolin) THEN
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     CALL qe_to_sirius_complex(v%ns_nc(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  ELSE
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat))
     CALL qe_to_sirius_real(v%ns(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  ENDIF
  ns_sirius(:,:,:,:) = 0.5 * ns_sirius(:,:,:,:)
  CALL sirius_set_hubbard_potential(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)

  DEALLOCATE(ns_sirius)
END SUBROUTINE qe_sirius_set_hubbard_potential

SUBROUTINE qe_sirius_get_hubbard_potential(v)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  IMPLICIT NONE

  COMPLEX(8), ALLOCATABLE :: ns_sirius(:, :, :, :)
  TYPE(scf_type), INTENT(out) :: v  ! the valence charge

  CALL sirius_get_hubbard_potential(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)

  IF (noncolin) THEN
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     ns_sirius(:,:,:,:) = 2.0 * ns_sirius(:,:,:,:)
     CALL sirius_to_qe_complex(ns_sirius(1, 1, 1, 1), v%ns_nc(1, 1, 1, 1))

  ELSE
     ALLOCATE(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat))
     ns_sirius(:,:,:,:) = 2.0 * ns_sirius(:,:,:,:)
     CALL sirius_to_qe_real(ns_sirius(1, 1, 1, 1), v%ns(1, 1, 1, 1))
  ENDIF

  DEALLOCATE(ns_sirius)
END SUBROUTINE qe_sirius_get_hubbard_potential


SUBROUTINE check_residuals(ik)
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE noncollin_module,     ONLY : noncolin, npol
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : okvan, nkb, vkb
  USE bp,                   ONLY : lelfield, bec_evcel
  USE becmod,               ONLY : bec_type, becp, calbec,&
                                   allocate_bec_type, deallocate_bec_type
  USE mp_bands,             ONLY : nproc_bgrp, intra_bgrp_comm, inter_bgrp_comm
IMPLICIT NONE
INTEGER, INTENT(in) :: ik
INTEGER npw, i, j, ig
COMPLEX(8), ALLOCATABLE :: hpsi(:,:), spsi(:,:), ovlp(:,:)
REAL(8) l2norm, l2norm_psi, max_err

ALLOCATE(hpsi(npwx*npol,nbnd))
ALLOCATE(spsi(npwx*npol,nbnd))
ALLOCATE(ovlp(nbnd, nbnd))

!  call allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
!do ik = 1, nks
!
!   if ( lsda ) current_spin = isk(ik)
   npw = ngk (ik)
!   !
!   if ( nks > 1 ) &
!      call get_buffer ( evc, nwordwfc, iunwfc, ik )
!
  !call g2_kin( ik )
!  !
!  ! ... more stuff needed by the hamiltonian: nonlocal projectors
!  !
!  if ( nkb > 0 ) call init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
!  call calbec(npw, vkb, evc, becp)

   CALL h_psi(npwx, npw, nbnd, evc, hpsi)
   IF (okvan) THEN
     CALL s_psi(npwx, npw, nbnd, evc, spsi)
   ELSE
     spsi = evc
   ENDIF

   DO j = 1, nbnd
     l2norm = 0
     l2norm_psi = 0
     DO ig = 1, npw
       l2norm = l2norm + ABS(hpsi(ig, j) - et(j, ik) * spsi(ig, j))**2
       l2norm_psi = l2norm_psi + REAL(CONJG(evc(ig, j)) * spsi(ig, j))
     ENDDO
     WRITE(*,*)'band: ', j, ', residual l2norm: ',SQRT(l2norm), ', psi norm: ',SQRT(l2norm_psi)

   ENDDO

   ovlp = 0
   max_err = 0
   DO i = 1, nbnd
     DO j = 1, nbnd
       DO ig = 1, npw
         ovlp(i, j) = ovlp(i, j) + CONJG(evc(ig, i)) * spsi(ig, j)
       ENDDO
       IF (i.EQ.j) ovlp(i, j) = ovlp(i, j) - 1.0
       IF (ABS(ovlp(i, j)).GT.max_err) THEN
         max_err = ABS(ovlp(i, j))
       ENDIF
       !if (abs(ovlp(i, j)).gt.1e-13) then
       !  write(*,*)'bands i, j:',i,j,', overlap: ', ovlp(i, j)
       !endif
     ENDDO
   ENDDO
   WRITE(*,*)'maximum overlap error: ',max_err


!enddo
!call deallocate_bec_type ( becp )

DEALLOCATE(hpsi, spsi, ovlp)

END SUBROUTINE check_residuals


END MODULE mod_sirius
