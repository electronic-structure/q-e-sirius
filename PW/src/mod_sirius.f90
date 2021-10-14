MODULE mod_sirius
USE input_parameters, ONLY : sirius_cfg, use_sirius_scf, use_sirius_nlcg
USE sirius
USE funct
USE mod_sirius_callbacks
USE mod_sirius_base
IMPLICIT NONE

CONTAINS

SUBROUTINE setup_sirius()
        USE constants, ONLY : RYTOEV
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
        & Hubbard_J0, U_projection, Hubbard_l, sc_at, at_sc, Hubbard_V, neighood, ldim_u
USE esm,                  ONLY : do_comp_esm
USE control_flags,        ONLY : iverbosity
USE Coul_cut_2D,          ONLY : do_cutoff_2D
USE mp_world, ONLY: mpime
IMPLICIT NONE
!
INTEGER :: dims(3), i, ia, iat, rank, ierr, ijv, li, lj, mb, nb, j, l,&
        ilast, ir, num_gvec, num_ranks_k, vt(3), iwf, num_kp, nmagd, atom_pair(2), &
        n_pair(2), l_pair(2), iwf1, iwf2, ia1, ia2, nt1, nt2, n1, n2, viz
REAL(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp, V_
REAL(8), ALLOCATABLE :: dion(:, :), qij(:,:,:), vloc(:), wk_tmp(:), xk_tmp(:,:)
INTEGER, ALLOCATABLE :: nk_loc(:)
INTEGER :: ih, jh, ijh, lmax_beta, nsymop
CHARACTER(LEN=1024) :: conf_str
INTEGER, EXTERNAL :: set_hubbard_l,set_hubbard_n
REAL(8), EXTERNAL :: hubbard_occ
REAL(8), PARAMETER :: spglib_tol=1e-4

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
           &"mixer" : {"beta" : ', F12.6, ', "max_history" : ', I4, ', "use_hartree" : false},&
           &"settings" : {"itsol_tol_scale" : [0.1, 0.95]},&
           &"control" : {"print_checksum" : false}}')
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
        &hubbard_correction=lda_plus_U, hubbard_correction_kind=lda_plus_u_kind, &
        &hubbard_full_orthogonalization=.true., &
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
 CALL sirius_set_callback_function(sctx, "ps_atomic_wf", C_FUNLOC(calc_atomic_wfc_radial_integrals))
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
! initialize atom types
IF (sirius_pwpp) THEN

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

      IF (ALLOCATED(upf(iat)%els)) then
         CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_atomic_wf", &
                 & upf(iat)%chi(1:msh(iat), iwf), msh(iat), l=l, occ=upf(iat)%oc(iwf), n=i, &
                 & orbital_label=upf(iat)%els(iwf))
      else
         CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_atomic_wf", &
                 & upf(iat)%chi(1:msh(iat), iwf), msh(iat), l=l, occ=upf(iat)%oc(iwf), n=i)
 endif
 ENDDO

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
  end do

  if (lda_plus_U) then


  DO iat = 1, nat
    !
    nt1 = ityp(iat)
    !
    IF (lda_plus_u_kind .EQ. 2) then
      IF (ldim_u(nt1).GT.0) THEN
        ! QE double count the links while sirius works with the links only. So
        ! we need a factor 1/2 for the coupling constants to take into account
        ! the double counting
        DO viz = 1, neighood(iat)%num_neigh
          atom_pair(1) = iat - 1
          ia2 = neighood(iat)%neigh(viz)
          atom_pair(2) = at_sc(ia2)%at - 1
          nt2 = ityp(atom_pair(2) + 1)
          !     sirius does not make any distinction
          ! between orbitals contributing to the U
          ! correction and orbitals with U = 0.
          ! Add them to the list of interacting
          ! orbitals if (V \neq 0) but (U = 0)
          ! terms contributing to the standard-standard term in QE language
          n_pair(1) = set_hubbard_n(upf(nt1)%psd)
          n_pair(2) = set_hubbard_n(upf(nt2)%psd)
          l_pair(1) = Hubbard_l(nt1)
          l_pair(2) = Hubbard_l(nt2)
          V_ = Hubbard_V(iat,ia2,1) * RYTOEV
          if (iat /= ia2) then
            call sirius_add_hubbard_atom_pair(sctx, &
                    atom_pair, &
                    at_sc(ia2)%n, &
                    n_pair, &
                    l_pair, &
                    V_)
            ! terms contributing to the standard-background term in QE language
            if (Abs(Hubbard_V(iat,ia2,2)) > 1e-8) then
              n2 = set_hubbard_n(upf(nt2)%psd)
              DO iwf = 1, upf(nt2)%nwfc
                if (n2 /= upf(nt2)%nchi(iwf)) then
                  n_pair(2) = upf(nt2)%nchi(iwf)
                  l_pair(2) = upf(nt2)%lchi(iwf)
                  V_ = Hubbard_V(iat,ia2,2)  * RYTOEV
                  call sirius_add_hubbard_atom_pair(sctx, &
                          atom_pair, &
                          at_sc(ia2)%n, &
                          n_pair, &
                          l_pair, &
                          V_)
                end if
              end DO
             end if
             ! terms contributing to the background-background  term in QE language
             if (Hubbard_V(iat,ia2,3) > 1e-8) then
               n1 = set_hubbard_n(upf(nt1)%psd)
               n2 = set_hubbard_n(upf(nt2)%psd)
               DO iwf1 = 1, upf(nt1)%nwfc
                 DO iwf2 = 1, upf(nt2)%nwfc
                   if ((n1 /= upf(nt1)%nchi(iwf)) .AND. (n2 /= upf(nt2)%nchi(iwf))) then
                     n_pair(1) = upf(nt1)%nchi(iwf)
                     n_pair(2) = upf(nt2)%nchi(iwf)
                     l_pair(1) = upf(nt1)%lchi(iwf)
                     l_pair(2) = upf(nt2)%lchi(iwf)
                     V_ = Hubbard_V(iat, ia2, 3)  * RYTOEV
                     call sirius_add_hubbard_atom_pair(sctx, &
                             atom_pair, &
                             at_sc(ia2)%n, &
                             n_pair, &
                             l_pair, &
                             V_)
                   end if
                 end DO
               end do
             END if
           else
             CALL sirius_set_atom_type_hubbard(sctx, &
                     & TRIM(atom_type(nt1)%label), &
                     & Hubbard_l(nt1), &
                     & set_hubbard_n(upf(nt1)%psd), &
                     & hubbard_occ ( upf(nt1)%psd ), &
                     & Hubbard_V(iat,ia2,1) * RYTOEV, &
                     & 0.0d0, &
                     & 0.0d0, &
                     & 0.0d0, &
                     & 0.0d0)
           end if
         END DO
       end if
     else
       IF (is_hubbard(iat)) THEN
         nt1 = ityp(iat)
         CALL sirius_set_atom_type_hubbard(sctx, &
                 & TRIM(atom_type(nt1)%label), &
                 & Hubbard_l(nt1), &
                 & set_hubbard_n(upf(nt1)%psd), &
                 & hubbard_occ ( upf(nt1)%psd ), &
                 & Hubbard_U(nt1) * RYTOEV, &
                 & Hubbard_J(1,nt1) * RYTOEV, &
                 & Hubbard_alpha(nt1) * RYTOEV, &
                 & Hubbard_beta(nt1) * RYTOEV, &
                 & Hubbard_J0(nt1) * RYTOEV)
       ENDIF
     END IF
   ENDDO ! iat
 endif

DO iat = 1, nsp
  ! update the hubbard information if needed otherwise do nothing
  CALL sirius_atom_type_update(sctx, TRIM(atom_type(iat)%label))
END DO
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
  !call sirius_reduce_coordinates(v1(1), v2(1), vt(1))
  v2 = v1
  IF (noncolin) THEN
          v1(1) = starting_magnetization(iat) * SIN(angle1(iat)) * COS(angle2(iat)) ! zv(iat) *
          v1(2) = starting_magnetization(iat) * SIN(angle1(iat)) * SIN(angle2(iat)) ! zv(iat) *
    v1(3) = starting_magnetization(iat) * COS(angle1(iat)) !zv(iat) *
  ELSE
    v1 = 0
    v1(3) = starting_magnetization(iat) * 3.0 !zv(iat) *
  ENDIF
  CALL sirius_add_atom(sctx, TRIM(atom_type(iat)%label), v2, v1)
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
  CALL sirius_dump_runtime_setup(sctx, "setup.json");
ENDIF

write(*,*)''
write(*,*)'=========================================='
write(*,*)'* SIRIUS simulation context finished     *'
write(*,*)'=========================================='


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
CALL sirius_stop_timer("setup_sirius")


write(*,*)''
write(*,*)'=============================================='
write(*,*)'* SIRIUS ground state structure initialized *'
write(*,*)'=============================================='

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
        CALL sirius_get_q_operator(sctx, TRIM(atom_type(iat)%label), ih, jh, ngm, mill, atom_type(iat)%qpw(:, ijh))
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
        CALL sirius_get_q_operator_matrix(sctx, TRIM(atom_type(iat)%label), qq_nt(1, 1, iat), nhm)
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

SUBROUTINE put_q_operator_matrix_to_sirius
USE uspp,       ONLY : qq_nt
USE uspp_param, ONLY : upf, nhm
USE ions_base,  ONLY : nsp, atm
IMPLICIT NONE
INTEGER iat

DO iat = 1, nsp
  IF (upf(iat)%tvanp) THEN
    CALL sirius_set_q_operator_matrix(sctx, TRIM(atom_type(iat)%label), qq_nt(1, 1, iat), nhm)
  ENDIF
ENDDO

END SUBROUTINE put_q_operator_matrix_to_sirius


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
INTEGER, ALLOCATABLE :: gvl(:,:)
INTEGER ig, ik, ik_, i, j, ispn, rank, ierr, nksmax
COMPLEX(8) z1
LOGICAL exst_file,exst_mem

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
