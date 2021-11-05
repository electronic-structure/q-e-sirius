MODULE mod_sirius
USE input_parameters, ONLY : sirius_cfg, use_sirius_scf, use_sirius_nlcg
USE sirius
USE funct
USE mod_sirius_callbacks
USE mod_sirius_base
IMPLICIT NONE

CONTAINS

SUBROUTINE setup_sirius()
USE cell_base, ONLY : alat, at, bg
USE ions_base, ONLY : tau, nsp, atm, zv, amass, ityp, nat
USE uspp_param, ONLY : upf, nhm, nh
USE atom, ONLY : rgrid, msh
USE fft_base, ONLY :  dfftp
USE klist, ONLY : nks, xk, nkstot, wk, degauss, ngauss
USE gvect, ONLY : ngm_g, ecutrho, ngm, mill
USE gvecw, ONLY : ecutwfc
USE control_flags, ONLY : gamma_only, diago_full_acc, mixing_beta, nmix
USE mp_pools, ONLY : inter_pool_comm, intra_pool_comm, npool
USE mp_images,        ONLY : nproc_image, intra_image_comm
USE mp, ONLY : mp_sum, mp_bcast
USE wvfct, ONLY : nbnd
USE parallel_include
USE sirius
USE input_parameters, ONLY : sirius_cfg, diago_david_ndim
USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
USE lsda_mod, ONLY : lsda, nspin, starting_magnetization
USE cell_base, ONLY : omega
USE symm_base,            ONLY : nosym, nsym
USE spin_orb,             ONLY : lspinorb
USE ldaU,                 ONLY : lda_plus_U, Hubbard_J, Hubbard_U, Hubbard_alpha, &
                               & Hubbard_beta, is_Hubbard, lda_plus_u_kind, &
                               & Hubbard_J0, U_projection, Hubbard_l
USE esm,                  ONLY : do_comp_esm
USE control_flags,        ONLY : iverbosity
USE Coul_cut_2D,          ONLY : do_cutoff_2D
USE mp_world, ONLY: mpime
USE klist,            ONLY : nkstot, nks, ngk, igk_k
IMPLICIT NONE
!
INTEGER :: dims(3), i, ia, iat, rank, ierr, ijv, li, lj, mb, nb, j, l,&
     ilast, ir, num_gvec, num_ranks_k, vt(3), iwf, num_kp, nmagd
REAL(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp
REAL(8), ALLOCATABLE :: dion(:, :), qij(:,:,:), vloc(:), wk_tmp(:), xk_tmp(:,:)
INTEGER, ALLOCATABLE :: nk_loc(:)
INTEGER :: ih, jh, ijh, lmax_beta, nsymop
CHARACTER(LEN=1024) :: conf_str
INTEGER, EXTERNAL :: set_hubbard_l,set_hubbard_n
REAL(8), EXTERNAL :: hubbard_occ
REAL(8), PARAMETER :: spglib_tol=1e-4
INTEGER, EXTERNAL :: global_kpoint_index

CALL sirius_start_timer("setup_sirius")

IF (do_cutoff_2D.OR.do_comp_esm) THEN
  use_veff_callback = .TRUE.
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
CALL sirius_create_context(intra_image_comm, sctx, fcomm_k=inter_pool_comm, fcomm_band=intra_pool_comm)
! create initial configuration dictionary in JSON
WRITE(conf_str, 10)diago_david_ndim, mixing_beta, nmix
10 FORMAT('{"parameters" : {"electronic_structure_method" : "pseudopotential", "use_scf_correction" : true}, &
           &"iterative_solver" : {"residual_tolerance" : 1e-6, "subspace_size" : ',I4,'}, &
           &"mixer" : {"beta" : ', F12.6, ', "max_history" : ', I4, ', "use_hartree" : true},&
           &"settings" : {"itsol_tol_scale" : [0.1, 0.95]}}')
! set initial parameters
CALL sirius_import_parameters(sctx, conf_str)
! set default verbosity
CALL sirius_set_parameters(sctx, verbosity=MIN(1, iverbosity))
! import config file
CALL sirius_import_parameters(sctx, TRIM(ADJUSTL(sirius_cfg)))

CALL sirius_get_parameters(sctx, electronic_structure_method=conf_str)
IF (TRIM(ADJUSTL(conf_str)) .EQ. "pseudopotential") THEN
  sirius_pwpp = .TRUE.
ELSE
  sirius_pwpp = .FALSE.
END IF

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
  &use_symmetry=(use_sirius_scf.OR.use_sirius_nlcg).AND..NOT.nosym, so_correction=lspinorb,&
  &pw_cutoff=SQRT(ecutrho), gk_cutoff=SQRT(ecutwfc),&
  &hubbard_correction=lda_plus_U, hubbard_correction_kind=lda_plus_u_kind,&
  &hubbard_orbitals=TRIM(ADJUSTL(U_projection)),fft_grid_size=dims, spglib_tol=spglib_tol)

! degauss is converted to Ha units
CALL sirius_set_parameters(sctx, smearing_width=degauss/2.d0)
SELECT CASE(ngauss)
  CASE (0)
    CALL sirius_set_parameters(sctx, smearing="gaussian")
  CASE(-1)
    CALL sirius_set_parameters(sctx, smearing="cold")
  CASE(-99)
    CALL sirius_set_parameters(sctx, smearing="fermi_dirac")
END SELECT

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
CALL sirius_set_callback_function(sctx, "ps_rho_ri", C_FUNLOC(calc_ps_rho_radial_integrals))
IF (use_veff_callback) THEN
  CALL sirius_set_callback_function(sctx, "veff", C_FUNLOC(calc_veff))
ENDIF
!CALL sirius_set_callback_function(sctx, "band_occ", C_FUNLOC(calc_band_occupancies))

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

write(100+rank,*)'nkstot=',nkstot, ' nks=',nks
do i = 1, nks
  write(100+rank,*)'ikloc=', i, ' ikglob=',global_kpoint_index(nkstot, i),' ngk=',ngk(i)
enddo
flush(100+rank)

IF (ALLOCATED(kpoint_index_map)) DEALLOCATE(kpoint_index_map)
ALLOCATE(kpoint_index_map(2, nkstot))
kpoint_index_map = 0
DO i = 1, nks
  kpoint_index_map(1, global_kpoint_index(nkstot, i)) = rank
  kpoint_index_map(2, global_kpoint_index(nkstot, i)) = i
END DO

call mpi_allreduce(MPI_IN_PLACE, kpoint_index_map, 2 * nkstot, MPI_INT, MPI_SUM, inter_pool_comm, ierr)
write(100+rank, *)'kpoint_index_map'
do i = 1, nkstot
  write(100+rank, *)i, kpoint_index_map(1, i), kpoint_index_map(2, i)
enddo
flush(100+rank)


IF (sirius_pwpp) THEN

  ! initialize atom types
  DO iat = 1, nsp

    ! add new atom type
     CALL sirius_add_atom_type(sctx, TRIM(atom_type(iat)%label), &
          & zn=NINT(zv(iat)+0.001d0), &
          & mass=amass(iat), &
          & spin_orbit=upf(iat)%has_so)

    ! set radial grid
    CALL sirius_set_atom_type_radial_grid(sctx, TRIM(atom_type(iat)%label), upf(iat)%mesh, upf(iat)%r)

    ! set beta-projectors
    DO i = 1, upf(iat)%nbeta
      l = upf(iat)%lll(i);
      IF (upf(iat)%has_so) THEN
        IF (upf(iat)%jjj(i) .LE. upf(iat)%lll(i)) THEN
          l = - upf(iat)%lll(i)
        ENDIF
      ENDIF
      CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "beta", &
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
      CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_atomic_wf", &
           & upf(iat)%chi(1:msh(iat), iwf), msh(iat), l=l, occ=upf(iat)%oc(iwf), n=i)
    ENDDO
  
    IF (is_hubbard(iat)) THEN
       CALL sirius_set_atom_type_hubbard(sctx, &
            & TRIM(atom_type(iat)%label), &
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
    CALL sirius_set_atom_type_dion(sctx, TRIM(atom_type(iat)%label), upf(iat)%nbeta, dion(1, 1))
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
            CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "q_aug",&
                                                     &upf(iat)%qfuncl(1:upf(iat)%kkbeta, ijv, l), upf(iat)%kkbeta,&
                                                     &l=l, idxrf1=i, idxrf2=j)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  
    IF (upf(iat)%tpawp) THEN
      DO i = 1, upf(iat)%nbeta
        CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ae_paw_wf",&
                                                 &upf(iat)%aewfc(1:upf(iat)%paw%iraug,i), upf(iat)%paw%iraug)
        CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_paw_wf",&
                                                 &upf(iat)%pswfc(1:upf(iat)%paw%iraug,i), upf(iat)%paw%iraug)
      ENDDO
      CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ae_paw_core",&
                                               &upf(iat)%paw%ae_rho_atc, upf(iat)%mesh)
  
      CALL sirius_set_atom_type_paw(sctx, TRIM(atom_type(iat)%label), upf(iat)%paw%core_energy / 2,&
                                   &upf(iat)%paw%oc, upf(iat)%nbeta)
    ENDIF
  
    ! set non-linear core correction
    IF (.TRUE.) THEN
      ALLOCATE(vloc(upf(iat)%mesh))
      vloc = 0.d0
      IF (ALLOCATED(upf(iat)%rho_atc)) THEN
        DO i = 1, msh(iat)
          vloc(i) = upf(iat)%rho_atc(i)
        ENDDO
      ENDIF
      CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_rho_core",&
                                               &vloc, upf(iat)%mesh)
      DEALLOCATE(vloc)
    ENDIF
  
    ! set total charge density of a free atom (to compute initial rho(r))
    CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_rho_total",&
                                             &upf(iat)%rho_at, upf(iat)%mesh)
  
    ! the hack is done in Modules/readpp.f90
    IF (.TRUE.) THEN
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
      CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "vloc", vloc, upf(iat)%mesh)
      DEALLOCATE(vloc)
    ENDIF
  ENDDO ! iat
ELSE
  ! initialize atom types with FP-LAPW species
  DO iat = 1, nsp
    ! add new atom type
     CALL sirius_add_atom_type(sctx, TRIM(atom_type(iat)%label), fname=TRIM(atom_type(iat)%label)//'.json')
  END DO
END IF

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
  IF (noncolin) THEN
    v2(1) = zv(iat) * starting_magnetization(iat) * SIN(angle1(iat)) * COS(angle2(iat))
    v2(2) = zv(iat) * starting_magnetization(iat) * SIN(angle1(iat)) * SIN(angle2(iat))
    v2(3) = zv(iat) * starting_magnetization(iat) * COS(angle1(iat))
  ELSE
    v2 = 0
    v2(3) = zv(iat) * starting_magnetization(iat)
  ENDIF
  CALL sirius_add_atom(sctx, TRIM(atom_type(iat)%label), v1, v2)
ENDDO

CALL put_xc_functional_to_sirius()

write(*,*)''
write(*,*)'=========================================='
write(*,*)'* initializing SIRIUS simulation context *'
write(*,*)'=========================================='

! initialize global variables/indices/arrays/etc. of the simulation
CALL sirius_initialize_context(sctx)

IF (sirius_pwpp) THEN
  DO iat = 1, nsp
    CALL sirius_get_num_beta_projectors(sctx, TRIM(atom_type(iat)%label), atom_type(iat)%num_beta_projectors)
  ENDDO
END IF

CALL sirius_get_parameters(sctx, num_sym_op=nsymop)
IF (nsymop .NE. nsym) THEN
  WRITE(*,*)
  WRITE(*,'("WARNING! Different number of symmetry operations: ", I4, " QE ", I4," spglib")')nsym, nsymop
  WRITE(*,*)
ENDIF

IF (mpime.eq.0) THEN
  CALL sirius_dump_runtime_setup(sctx, "setup.json")
ENDIF



! get number of g-vectors of the dense fft grid
CALL sirius_get_num_gvec(sctx, num_gvec)

IF (.NOT.((num_gvec .EQ. ngm_g) .OR. (num_gvec * 2 - 1 .EQ. ngm_g))) THEN
  WRITE(*,*)
  WRITE(*,'("Error: wrong number of G-vectors; QE: ", I6, ", SIRIUS: ", I6)')ngm_g, num_gvec
  WRITE(*,*)
  STOP
END IF

! create k-point set
! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice and
!          without x2 multiplication for the lsda case
CALL sirius_create_kset(sctx, num_kpoints, kpoints, wkpoints, .FALSE., ks_handler)

! create ground-state class
CALL sirius_create_ground_state(ks_handler, gs_handler)

CALL sirius_stop_timer("setup_sirius")

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
  DO ik = 1, nkstot
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
  DO ik = 1, nkstot
    wg(:, ik) = bnd_occ(:, global_kpoint_index(nkstot, ik)) / maxocc * wk(ik)
  ENDDO

  DEALLOCATE(bnd_occ)

END SUBROUTINE get_band_occupancies_from_sirius


SUBROUTINE get_wave_functions_from_sirius
USE klist,            ONLY : nkstot, nks, ngk, igk_k
USE gvect,            ONLY : mill
USE buffers,          ONLY : save_buffer, open_buffer
USE io_files,         ONLY : iunwfc, nwordwfc
USE control_flags,    ONLY : io_level
USE bp,               ONLY : lelfield
USE noncollin_module, ONLY : npol
USE wvfct,            ONLY : npwx, nbnd
USE wavefunctions,    ONLY : evc
USE lsda_mod,         ONLY : isk, lsda
USE mp_pools,         ONLY : inter_pool_comm
USE parallel_include
IMPLICIT NONE
INTEGER, EXTERNAL :: global_kpoint_index
INTEGER, ALLOCATABLE :: vgl(:,:)
INTEGER ig, ik, ik_, ik1, i, j, ispn, rank, ierr, nksmax, ikloc
COMPLEX(8) z1
LOGICAL exst_file,exst_mem

! rank of communicator that distributes k-points
CALL mpi_comm_rank(inter_pool_comm, rank, ierr)

ALLOCATE(vgl(3, npwx))

DO ik = 1, nkstot
  ik1 = MOD(ik - 1, num_kpoints) + 1
  ispn = 1
  IF (ik .GT. num_kpoints) THEN
    ispn = 2
  ENDIF
  IF (kpoint_index_map(1, ik) .EQ. rank) THEN
    ikloc = kpoint_index_map(2, ik)
    DO ig = 1, ngk(ikloc)
      vgl(:,ig) = mill(:, igk_k(ig, ikloc))
    ENDDO
    CALL sirius_get_wave_functions_v2( ks_handler, vkl=kpoints(:, ik1), spin=ispn, num_gvec_loc=ngk(ikloc), &
                                     & gvec_loc=vgl(1, 1), evec=evc(1, 1), ld=npwx, num_spin_comp=npol )
    IF (nks > 1 .OR. lelfield) THEN
      CALL save_buffer ( evc, nwordwfc, iunwfc, ikloc )
    ENDIF
  ELSE
    CALL sirius_get_wave_functions_v2( ks_handler )
  ENDIF

  CALL mpi_barrier(inter_pool_comm, ierr)

ENDDO

!CALL mpi_allreduce(nks, nksmax, 1, MPI_INTEGER, MPI_MAX, inter_pool_comm, ierr)
!
!!CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
!
!DO ik = 1, nksmax
!  IF (ik.LE.nks) THEN
!    DO ig = 1, ngk(ik)
!      gvl(:,ig) = mill(:, igk_k(ig, ik))
!    ENDDO
!    !
!    ik_ = global_kpoint_index(nkstot, ik)
!    ispn = isk(ik)
!    IF (lsda.AND.ispn.EQ.2) THEN
!      ik_ = ik_ - nkstot / 2
!    ENDIF
!    CALL sirius_get_wave_functions(ks_handler, ik_, ispn, ngk(ik), gvl(1, 1), evc(1, 1), npwx, npol)
!    IF (nks > 1 .OR. lelfield) THEN
!      CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
!    ENDIF
!  ELSE
!    CALL sirius_get_wave_functions(ks_handler, -1, -1, -1, -1, z1, -1, -1)
!  ENDIF
!ENDDO

DEALLOCATE(vgl)

END SUBROUTINE get_wave_functions_from_sirius


SUBROUTINE put_xc_functional_to_sirius
USE xc_lib
IMPLICIT NONE
INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac

  iexch  = xclib_get_id( 'LDA', 'EXCH' )
  icorr  = xclib_get_id( 'LDA', 'CORR' )
  igcx   = xclib_get_id( 'GGA', 'EXCH' )
  igcc   = xclib_get_id( 'GGA', 'CORR' )
  imeta  = xclib_get_id( 'MGGA','EXCH' )
  imetac = xclib_get_id( 'MGGA','CORR' )

  IF (imeta.NE.0.OR.imetac.NE.0) THEN
   STOP ("interface for meta-XC functional is not implemented")
  ENDIF

  IF (iexch.NE.0.AND.igcx.EQ.0) THEN
   SELECT CASE(iexch)
   CASE(0)
   CASE(1)
     CALL sirius_add_xc_functional(sctx, "XC_LDA_X")
   CASE default
     STOP ("interface for this exchange functional is not implemented")
   END SELECT
  ENDIF

  IF (iexch.NE.0.AND.igcx.NE.0) THEN
   SELECT CASE(igcx)
   CASE(0)
   CASE(2)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PW91")
   CASE(3)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE")
   CASE(10)
     CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE_SOL")
   CASE default
     WRITE(*,*)igcx
     STOP ("interface for this gradient exchange functional is not implemented")
   END SELECT
  ENDIF

  IF (icorr.NE.0.AND.igcc.EQ.0) THEN
   SELECT CASE(icorr)
   CASE(0)
   CASE(1)
     CALL sirius_add_xc_functional(sctx, "XC_LDA_C_PZ")
   CASE(4)
     CALL sirius_add_xc_functional(sctx, "XC_LDA_C_PW")
   CASE default
     STOP ("interface for this correlation functional is not implemented")
   END SELECT
  ENDIF

  IF (icorr.NE.0.AND.igcc.NE.0) THEN
   SELECT CASE(igcc)
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


END MODULE mod_sirius
