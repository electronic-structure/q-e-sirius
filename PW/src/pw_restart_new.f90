!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE pw_restart_new
!----------------------------------------------------------------------------
  !
  ! ... New PWscf I/O using xml schema and hdf5 binaries
  ! ... Parallel execution: the xml file is written by one processor only
  ! ... ("ionode_id"), read by all processors ;
  ! ... the wavefunction files are written / read by one processor per pool,
  ! ... collected on / distributed to all other processors in pool
  !
  USE KINDS,        ONLY: DP
  USE qes_types_module
  !USE qes_libs_module, ONLY: qes_write, qes_reset, qes_init
  USE qes_write_module, ONLY: qes_write
  USE qes_reset_module, ONLY: qes_reset 
  USE qes_init_module, ONLY: qes_init
  USE qexsd_module, ONLY: qexsd_init_schema, qexsd_openschema, qexsd_closeschema,      &
                          qexsd_init_convergence_info, qexsd_init_algorithmic_info,    & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure,      &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure,          &
                          qexsd_init_dipole_info, qexsd_init_total_energy,             &
                          qexsd_init_forces,qexsd_init_stress, qexsd_xf,               &
                          qexsd_init_outputElectricField,                              &
                          qexsd_input_obj, qexsd_occ_obj, qexsd_smear_obj,             &
                          qexsd_init_outputPBC, qexsd_init_gate_info, qexsd_init_hybrid,&
                          qexsd_init_dftU, qexsd_init_vdw
  USE io_global, ONLY : ionode, ionode_id
  USE io_files,  ONLY : iunpun, xmlpun_schema, prefix, tmp_dir, postfix
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  PRIVATE
  PUBLIC :: pw_write_schema, pw_write_binaries, &
       pw_readschema_file, init_vars_from_schema, read_collected_to_evc
  !
  CONTAINS
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_schema( what, wf_collect )
      !------------------------------------------------------------------------
      !
      ! what = 'init-config': write only variables that are known after the 
      !                       initial steps of initialization (e.g. structure)
      ! wf_collect = T  if final wavefunctions in portable format are written,
      !              F  if wavefunctions are either not written or are written
      !                 in binary non-portable form (for checkpointing)
      !                 NB: wavefunctions are not written here in any case
      !
      USE control_flags,        ONLY : istep, conv_ions, &
                                       lscf, gamma_only, &
                                       tqr, tq_smoothing, tbeta_smoothing, &
                                       noinv, smallmem, &
                                       llondon, lxdm, ts_vdw, scf_error, n_scf_steps
      USE constants,            ONLY : e2
      USE realus,               ONLY : real_space
      USE uspp,                 ONLY : okvan
      USE paw_variables,        ONLY : okpaw
      USE uspp_param,           ONLY : upf
      USE global_version,       ONLY : version_number
      USE cell_base,            ONLY : at, bg, alat, ibrav
      USE gvect,                ONLY : ig_l2g
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, zv
      USE noncollin_module,     ONLY : noncolin, npol
      USE io_files,             ONLY : nwordwfc, iunwfc, psfile
      USE buffers,              ONLY : get_buffer
      USE wavefunctions, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, wk, &
                                       lgauss, ngauss, smearing, degauss, nelec, &
                                       two_fermi_energies, nelup, neldw, tot_charge
      USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE gvect,                ONLY : ngm, ngm_g, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,                ONLY : ngms_g, dual
      USE fft_base,             ONLY : dffts
      USE wvfct,                ONLY : npwx, et, wg, nbnd
      USE ener,                 ONLY : ef, ef_up, ef_dw, vtxc, etxc, ewld, etot, &
                                       ehart, eband, demet, edftd3, elondon, exdm
      USE tsvdw_module,         ONLY : EtsvdW
      USE gvecw,                ONLY : ecutwfc
      USE fixed_occ,            ONLY : tfixed_occ, f_inp
      USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, U_projection, &
                                       Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, &
                                       Hubbard_alpha, Hubbard_J0, Hubbard_beta,&
                                       is_hubbard
      USE spin_orb,             ONLY : lspinorb, domag
      USE symm_base,            ONLY : nrot, nsym, invsym, s, ft, irt, &
                                       t_rev, sname, time_reversal, no_t_rev,&
                                       spacegroup
      USE lsda_mod,             ONLY : nspin, isk, lsda, starting_magnetization, magtot, absmag
      USE noncollin_module,     ONLY : angle1, angle2, i_cons, mcons, bfield, magtot_nc, &
                                       lambda
      USE ions_base,            ONLY : amass
      USE funct,                ONLY : get_dft_short, get_inlc, get_nonlocc_name, dft_is_nonlocc
      USE scf,                  ONLY : rho
      USE force_mod,            ONLY : lforce, sumfor, force, sigma, lstres
      USE extfield,             ONLY : tefield, dipfield, edir, etotefield, &
                                       emaxpos, eopreg, eamp, el_dipole, ion_dipole,&
                                       gate, zgate, relaxz, block, block_1,&
                                       block_2, block_height, etotgatefield ! TB
      USE mp,                   ONLY : mp_sum
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE funct,                ONLY : get_exx_fraction, dft_is_hybrid, &
                                       get_gau_parameter, &
                                       get_screening_parameter, exx_is_active
      USE exx_base,             ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut
      USE exx,                  ONLY : ecutfock, local_thr 
      USE london_module,        ONLY : scal6, lon_rcut, c6_i
      USE xdm_module,           ONLY : xdm_a1=>a1i, xdm_a2=>a2i
      USE tsvdw_module,         ONLY : vdw_isolated, vdw_econv_thr
      USE input_parameters,     ONLY : verbosity, calculation, ion_dynamics, starting_ns_eigenvalue, &
                                       vdw_corr, london, k_points, assume_isolated, &  
                                       input_parameters_occupations => occupations, dftd3_threebody, &
                                       dftd3_version
      USE bp,                   ONLY : lelfield, lberry, el_pol, ion_pol
      !
      USE rap_point_group,      ONLY : elem, nelem, name_class
      USE rap_point_group_so,   ONLY : elem_so, nelem_so, name_class_so
      USE bfgs_module,          ONLY : bfgs_get_n_iter
      USE qexsd_module,         ONLY : qexsd_bp_obj, qexsd_start_k_obj
      USE qexsd_input,          ONLY : qexsd_init_k_points_ibz, &
              qexsd_init_occupations, qexsd_init_smearing
      USE fcp_variables,        ONLY : lfcpopt, lfcpdyn, fcp_mu  
      USE io_files,             ONLY : pseudo_dir
      USE control_flags,        ONLY : conv_elec, conv_ions, ldftd3
      USE input_parameters,     ONLY :  ts_vdw_econv_thr, ts_vdw_isolated
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      LOGICAL, INTENT(IN) :: wf_collect
      !
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname
      INTEGER               :: i, ig, ngg, ipol
      INTEGER               :: npwx_g, ispin, inlc
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      LOGICAL               :: lwfc, lrho, lxsd, occupations_are_fixed
      INTEGER                  :: iclass, isym, ielem
      CHARACTER(LEN=15)        :: symop_2_class(48)
      LOGICAL                  :: opt_conv_ispresent, dft_is_vdw, empirical_vdw
      INTEGER                  :: n_opt_steps, n_scf_steps_, h_band
      REAL(DP)                 :: h_energy
      TYPE(gateInfo_type),TARGET      :: gate_info_temp
      TYPE(gateInfo_type),POINTER     :: gate_info_ptr => NULL()
      TYPE(dipoleOutput_type),TARGET  :: dipol_obj 
      TYPE(dipoleOutput_type),POINTER :: dipol_ptr  => NULL()
      TYPE(BerryPhaseOutput_type),  POINTER :: bp_obj_ptr => NULL()
      TYPE(hybrid_type), POINTER            :: hybrid_obj => NULL()
      TYPE(vdW_type), POINTER               :: vdw_obj => NULL()
      TYPE(dftU_type), POINTER              :: dftU_obj => NULL() 
      REAL(DP), TARGET                      :: lumo_tmp, ef_targ, dispersion_energy_term 
      REAL(DP), POINTER                     :: lumo_energy => NULL(), ef_point => NULL()
      REAL(DP), ALLOCATABLE                 :: ef_updw(:)
      !
      !
      !
      TYPE(output_type) :: output
      REAL(DP),POINTER    :: degauss_, demet_, efield_corr, potstat_corr,  gatefield_corr  
      LOGICAL, POINTER    :: optimization_has_converged => NULL() 
      LOGICAL, TARGET     :: conv_opt  
      LOGICAL             :: scf_has_converged 
      INTEGER             :: itemp = 1
      REAL(DP),ALLOCATABLE :: london_c6_(:), bp_el_pol(:), bp_ion_pol(:), U_opt(:), J0_opt(:), alpha_opt(:), &
                              J_opt(:,:), beta_opt(:) 
      CHARACTER(LEN=3),ALLOCATABLE :: species_(:)
      CHARACTER(LEN=20),TARGET   :: dft_nonlocc_
      INTEGER,TARGET             :: dftd3_version_
      CHARACTER(LEN=20),TARGET   :: vdw_corr_
      CHARACTER(LEN=20),POINTER  :: non_local_term_pt =>NULL(), vdw_corr_pt=>NULL()
      REAL(DP),TARGET            :: temp(20), lond_rcut_, lond_s6_, ts_vdw_econv_thr_, xdm_a1_, xdm_a2_, ectuvcut_,&
                                    scr_par_, loc_thr_  
      REAL(DP),POINTER           :: vdw_term_pt =>NULL(), ts_thr_pt=>NULL(), london_s6_pt=>NULL(),&
                                    london_rcut_pt=>NULL(), xdm_a1_pt=>NULL(), xdm_a2_pt=>NULL(), &
                                    ts_vdw_econv_thr_pt=>NULL(), ectuvcut_opt=>NULL(), scr_par_opt=>NULL(), &
                                    loc_thr_p => NULL() 
      LOGICAL,TARGET             :: dftd3_threebody_, ts_vdw_isolated_
      LOGICAL,POINTER            :: ts_isol_pt=>NULL(), dftd3_threebody_pt=>NULL(), ts_vdw_isolated_pt =>NULL()
      INTEGER,POINTER            :: dftd3_version_pt => NULL() 

      NULLIFY( degauss_, demet_, efield_corr, potstat_corr, gatefield_corr )

      !
      ! Global PW dimensions need to be properly computed, reducing across MPI tasks
      ! If local PW dimensions are not available, set to 0
      !
      ALLOCATE( ngk_g( nkstot ) )
      ngk_g(:) = 0
      IF ( ALLOCATED (ngk) ) THEN
         ngk_g(1:nks) = ngk(:)
         CALL mp_sum( ngk_g(1:nks), intra_bgrp_comm )
         CALL ipoolrecover( ngk_g, 1, nkstot, nks )
      END IF
      ! BEWARE: only the first pool has ngk_g for all k-points
      !
      ! ... compute the maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      ! 
      ! 
      ! XML descriptor
      ! 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix
      !
      CALL qexsd_init_schema( iunpun )
      !
      !
      IF ( ionode ) THEN  
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL qexsd_openschema(TRIM( dirname ) // TRIM( xmlpun_schema ), 'PWSCF' )
         output%tagname="output"
         output%lwrite = .TRUE.
         output%lread  = .TRUE.
         !
!-------------------------------------------------------------------------------
! ... CONVERGENCE_INFO
!-------------------------------------------------------------------------------
         SELECT CASE (TRIM( calculation )) 
            CASE ( "relax","vc-relax" )
                conv_opt = conv_ions  
                optimization_has_converged  => conv_opt
                IF (TRIM( ion_dynamics) == 'bfgs' ) THEN 
                    n_opt_steps = bfgs_get_n_iter('bfgs_iter ') 
                ELSE 
                    n_opt_steps = istep 
                END IF 
                scf_has_converged = conv_elec 
                n_scf_steps_ = n_scf_steps
            CASE ("nscf", "bands" )
                n_opt_steps = 0
                scf_has_converged = .FALSE. 
                n_scf_steps_ = 1
            CASE default
                n_opt_steps        = 0 
                n_scf_steps_ = n_scf_steps
         END SELECT
         ! 
            call qexsd_init_convergence_info(output%convergence_info,   &
                        SCf_HAS_CONVERGED = scf_has_converged, &
                        OPTIMIZATION_HAS_CONVERGED = optimization_has_converged,& 
                        N_SCF_STEPS = n_scf_steps_, SCF_ERROR=scf_error,&
                        N_OPT_STEPS = n_opt_steps, GRAD_NORM = sumfor)
            output%convergence_info_ispresent = .TRUE.
         !
            
!-------------------------------------------------------------------------------
! ... ALGORITHMIC_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_algorithmic_info(output%algorithmic_info, &
              REAL_SPACE_BETA = real_space, REAL_SPACE_Q=tqr , USPP=okvan, PAW=okpaw)
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_SPECIES
!-------------------------------------------------------------------------------
         !
         ! while amass's are always present, starting_mag should not be passed
         ! for nspin==1 or contrained magnetization calculations
         !
         IF (noncolin) THEN
            CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm, psfile, &
                 amass, STARTING_MAGNETIZATION = starting_magnetization, &
                 ANGLE1=angle1, ANGLE2=angle2)
         ELSE IF (nspin==2) THEN 
            CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm, psfile, &
                 amass, STARTING_MAGNETIZATION=starting_magnetization)
         ELSE 
            CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm,psfile, &
                 amass)
         END IF
         output%atomic_species%pseudo_dir = TRIM(pseudo_dir)
         output%atomic_species%pseudo_dir_ispresent = .TRUE.
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_STRUCTURE
!-------------------------------------------------------------------------------
         !         
         CALL qexsd_init_atomic_structure(output%atomic_structure, nsp, atm, ityp, &
              nat, alat*tau, alat, alat*at(:,1), alat*at(:,2), alat*at(:,3), ibrav)
         !
!-------------------------------------------------------------------------------
! ... SYMMETRIES
!-------------------------------------------------------------------------------
         !
         symop_2_class="not found"
         IF (TRIM (verbosity) == 'medium' .OR. TRIM(verbosity) == 'high') THEN
            IF ( noncolin )  THEN 
               symmetries_so_loop:DO isym = 1, nrot 
                  classes_so_loop:DO iclass = 1, 24
                     elements_so_loop:DO ielem=1, nelem_so(iclass)
                        IF ( elem_so(ielem,iclass) == isym) THEN 
                           symop_2_class(isym) = name_class_so(iclass)
                           EXIT symmetries_so_loop
                        END IF
                     END DO elements_so_loop 
                     END DO classes_so_loop
               END DO symmetries_so_loop
            !
            ELSE
               symmetries_loop:DO isym = 1, nrot
                  classes_loop:DO iclass = 1, 12
                     elements_loop:DO ielem=1, nelem (iclass)
                        IF ( elem(ielem,iclass) == isym) THEN
                           symop_2_class(isym) = name_class(iclass)
                           EXIT classes_loop
                        END IF
                     END DO elements_loop
                  END DO classes_loop
               END DO symmetries_loop
            END IF
         END IF
         CALL qexsd_init_symmetries(output%symmetries, nsym, nrot, spacegroup,&
              s, ft, sname, t_rev, nat, irt,symop_2_class(1:nrot), verbosity, &
              noncolin)
         output%symmetries_ispresent=.TRUE. 
         !
!-------------------------------------------------------------------------------
! ... BASIS SET
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_basis_set(output%basis_set, gamma_only, ecutwfc/e2, ecutwfc*dual/e2, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, dffts%nr1, dffts%nr2, dffts%nr3, &
              .FALSE., dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, ngms_g, npwx_g, &
              bg(:,1), bg(:,2), bg(:,3) )
         !
!-------------------------------------------------------------------------------
! ... DFT
!-------------------------------------------------------------------------------
         !
         IF (dft_is_hybrid() ) THEN 
            ALLOCATE ( hybrid_obj)
            IF (get_screening_parameter() > 0.0_DP) THEN
               scr_par_ = get_screening_parameter() 
               scr_par_opt=> scr_par_ 
            END IF 
            IF (ecutvcut > 0.0_DP) THEN 
               ectuvcut_ = ecutvcut/e2 
               ectuvcut_opt => ectuvcut_
            END IF 
            IF ( local_thr > 0._DP) THEN 
               loc_thr_ = local_thr 
               loc_thr_p => loc_thr_ 
            END IF 
            CALL qexsd_init_hybrid(hybrid_obj, DFT_IS_HYBRID = .TRUE., NQ1 = nq1 , NQ2 = nq2, NQ3 =nq3, ECUTFOCK = ecutfock/e2, &
                                   EXX_FRACTION = get_exx_fraction(), SCREENING_PARAMETER = scr_par_opt, &
                                   EXXDIV_TREATMENT = exxdiv_treatment, X_GAMMA_EXTRAPOLATION = x_gamma_extrapolation,&
                                   ECUTVCUT = ectuvcut_opt, LOCAL_THR = loc_thr_p )
         END IF 

         empirical_vdw = (llondon .OR. ldftd3 .OR. lxdm .OR. ts_vdw )
         dft_is_vdw = dft_is_nonlocc() 
         IF ( dft_is_vdw .OR. empirical_vdw ) THEN 
            ALLOCATE (vdw_obj)
            IF ( empirical_vdw) THEN
                vdw_term_pt => dispersion_energy_term
                vdw_corr_ = TRIM(vdw_corr)
                vdw_corr_pt => vdw_corr_
                IF (llondon ) THEN
                    dispersion_energy_term = elondon/e2
                    lond_s6_ = scal6
                    london_s6_pt => lond_s6_
                    lond_rcut_ = lon_rcut
                    london_rcut_pt => lond_rcut_
                    IF (ANY( c6_i(1:nsp) .NE. -1._DP )) THEN
                       ALLOCATE (london_c6_(nsp), species_(nsp))
                       london_c6_(1:nsp) = c6_i(1:nsp)
                       species_(1:nsp)  = atm(1:nsp)
                   END IF
                   !
                ELSE IF ( lxdm ) THEN
                    dispersion_energy_term = exdm/e2
                    xdm_a1_ = xdm_a1
                    xdm_a1_pt => xdm_a1_
                    xdm_a2_ = xdm_a2
                    xdm_a2_pt => xdm_a2_
                    !
                ELSE IF ( ldftd3) THEN
                    dispersion_energy_term = edftd3/e2
                    dftd3_version_ = dftd3_version
                    dftd3_version_pt => dftd3_version_
                    dftd3_threebody_ = dftd3_threebody
                    dftd3_threebody_pt => dftd3_threebody_
                ELSE IF ( ts_vdw ) THEN
                    dispersion_energy_term = 2._DP * EtsvdW/e2
                    ts_vdw_isolated_ = ts_vdw_isolated
                    ts_vdw_isolated_pt => ts_vdw_isolated_
                    ts_vdw_econv_thr_ = ts_vdw_econv_thr
                    ts_vdw_econv_thr_pt => ts_vdw_econv_thr_
                END IF
            END IF 
            IF (dft_is_vdw) THEN
                dft_nonlocc_ = TRIM(get_nonlocc_name())
                non_local_term_pt => dft_nonlocc_
            END IF
            CALL qexsd_init_vdw(vdw_obj, non_local_term_pt, vdw_corr_pt, vdw_term_pt, &
                                ts_thr_pt, ts_isol_pt, london_s6_pt, LONDON_C6 = london_c6_, &
                                LONDON_RCUT =   london_rcut_pt, XDM_A1 = xdm_a1_pt, XDM_A2 = xdm_a2_pt,&
                                 DFTD3_VERSION = dftd3_version_pt, DFTD3_THREEBODY = dftd3_threebody_pt)
         END IF 
         IF ( lda_plus_u) THEN 
            ALLOCATE (dftU_obj)  
            CALL check_and_allocate(U_opt, Hubbard_U)
            CALL check_and_allocate(J0_opt, Hubbard_J0) 
            CALL check_and_allocate(alpha_opt, Hubbard_alpha) 
            CALL check_and_allocate(beta_opt, Hubbard_beta) 
            IF ( ANY(Hubbard_J(:,1:nsp) /= 0.0_DP)) THEN
               ALLOCATE (J_opt(3,nsp)) 
               J_opt(:, 1:nsp) = Hubbard_J(:, 1:nsp) 
            END IF 
            CALL qexsd_init_dftU (dftU_obj, NSP = nsp, SPECIES = atm(1:nsp), ITYP = ityp(1:nat),                     &
                                  IS_HUBBARD = is_hubbard, PSD = upf(1:nsp)%psd, NONCOLIN = noncolin, U =U_opt, &
                                  LDA_PLUS_U_KIND = lda_plus_u_kind, U_PROJECTION_TYPE = U_projection,               &
                                  J0 = J0_opt, alpha = alpha_opt, beta = beta_opt, J = J_opt,        & 
                                  starting_ns = starting_ns_eigenvalue, Hub_ns = rho%ns, Hub_ns_nc = rho%ns_nc ) 
         END IF 
         dft_name = get_dft_short()
         inlc = get_inlc()
         !
         CALL qexsd_init_dft  (output%dft, dft_name, hybrid_obj, vdw_obj, dftU_obj)
         IF (ASSOCIATED (hybrid_obj)) THEN
            CALL qes_reset(hybrid_obj) 
            DEALLOCATE (hybrid_obj) 
         END IF 
         IF (ASSOCIATED (vdw_obj)) THEN
            CALL qes_reset(vdw_obj) 
            DEALLOCATE (vdw_obj) 
         END IF 
         IF (ASSOCIATED (dftU_obj)) THEN 
            CALL qes_reset( dftU_obj) 
            DEALLOCATE (dftU_obj) 
         END IF 
         !
!-------------------------------------------------------------------------------
! ... PERIODIC BOUNDARY CONDITIONS 
!-------------------------------------------------------------------------------
         !
         IF (TRIM( assume_isolated ) .EQ. "2D" ) THEN
            output%boundary_conditions_ispresent=.TRUE.
            CALL  qexsd_init_outputPBC(output%boundary_conditions, assume_isolated)
          ENDIF
         !
!-------------------------------------------------------------------------------
! ... MAGNETIZATION
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_magnetization(output%magnetization, lsda, noncolin, lspinorb, &
              magtot, magtot_nc, absmag, domag )
         !

!--------------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------------
         !
         ! skip if not yet computed
         !
         IF ( TRIM(what) == "init-config" ) GO TO 10
         !
         IF (TRIM(input_parameters_occupations) == 'fixed') THEN 
            occupations_are_fixed = .TRUE.
            CALL get_homo_lumo( h_energy, lumo_tmp)
            h_energy = h_energy/e2
            IF ( lumo_tmp .LT. 1.d+6 ) THEN
                lumo_tmp = lumo_tmp/e2
                lumo_energy => lumo_tmp
            END IF
         ELSE 
            occupations_are_fixed = .FALSE. 
         END IF
         IF (nks_start == 0 .AND. nk1*nk2*nk3 > 0 ) THEN 
            CALL qexsd_init_k_points_ibz(qexsd_start_k_obj, "automatic", calculation, &
                 nk1, nk2, nk3, k1, k2, k3, nks_start, xk_start, wk_start, alat, at(:,1), .TRUE.)
         ELSE
            CALL qexsd_init_k_points_ibz(qexsd_start_k_obj, k_points, calculation, &
                                nk1, nk2, nk3, k1, k2, k3, nks_start, xk_start, wk_start, alat, at(:,1), .TRUE.)
         END IF
         qexsd_start_k_obj%tagname = 'starting_kpoints'
         IF ( TRIM (qexsd_input_obj%tagname) == 'input') THEN 
            qexsd_occ_obj = qexsd_input_obj%bands%occupations
         ELSE 
            CALL qexsd_init_occupations ( qexsd_occ_obj, input_parameters_occupations, nspin)
         END IF 
         qexsd_occ_obj%tagname = 'occupations_kind' 
         IF ( two_fermi_energies ) THEN
            ALLOCATE ( ef_updw (2) )
               IF (TRIM(input_parameters_occupations) == 'fixed') THEN  
                  ef_updw(1)  = MAXVAL(et(INT(nelup),1:nkstot/2))/e2
                  ef_updw (2)  = MAXVAL(et(INT(neldw),nkstot/2+1:nkstot))/e2 
               ELSE 
                  ef_updw = [ef_up/e2, ef_dw/e2]
               END IF
         ELSE
                ef_targ = ef/e2
                ef_point => ef_targ
         END IF


         IF (TRIM(input_parameters_occupations) == 'smearing' ) THEN
            IF (TRIM(qexsd_input_obj%tagname) == 'input') THEN 
               qexsd_smear_obj = qexsd_input_obj%bands%smearing
            ELSE 
               CALL qexsd_init_smearing(qexsd_smear_obj, smearing, degauss)
            END IF  
            !  
            CALL qexsd_init_band_structure(  output%band_structure,lsda,noncolin,lspinorb, nelec, natomwfc, &
                                 et, wg, nkstot, xk, ngk_g, wk, SMEARING = qexsd_smear_obj,  &
                                 STARTING_KPOINTS = qexsd_start_k_obj, OCCUPATIONS_KIND = qexsd_occ_obj, &
                                 WF_COLLECTED = wf_collect, NBND = nbnd, FERMI_ENERGY = ef_point, EF_UPDW = ef_updw )
            CALL qes_reset (qexsd_smear_obj)
         ELSE     
            CALL  qexsd_init_band_structure(output%band_structure,lsda, noncolin,lspinorb, nelec, natomwfc, &
                                et, wg, nkstot, xk, ngk_g, wk, &
                                STARTING_KPOINTS =  qexsd_start_k_obj, OCCUPATIONS_KIND = qexsd_occ_obj,&
                                WF_COLLECTED = wf_collect , NBND = nbnd, HOMO = h_energy, LUMO = lumo_energy , &
                                EF_UPDW = ef_updw )
         END IF 
         CALL qes_reset (qexsd_start_k_obj)
         CALL qes_reset (qexsd_occ_obj)
         !
!-------------------------------------------------------------------------------------------
! ... TOTAL ENERGY
!-------------------------------------------------------------------------------------------
         !
         IF ( degauss > 0.0d0 ) THEN
            !
            itemp = itemp + 1 
            temp(itemp)  = degauss/e2
            degauss_ => temp(itemp)
            !
            itemp = itemp+1
            temp(itemp)   = demet/e2
            demet_ => temp(itemp) 
         END IF
         IF ( tefield ) THEN 
            itemp = itemp+1 
            temp(itemp) = etotefield/e2
            efield_corr => temp(itemp) 
         END IF
         IF (lfcpopt .OR. lfcpdyn ) THEN 
            itemp = itemp +1 
            temp(itemp) = ef * tot_charge/e2
            potstat_corr => temp(itemp) 
            output%FCP_tot_charge_ispresent = .TRUE.
            output%FCP_tot_charge = tot_charge
            output%FCP_force_ispresent = .TRUE.
            !FIXME ( decide what units to use here ) 
            output%FCP_force = fcp_mu - ef 
         END IF 
         IF ( gate) THEN
            itemp = itemp + 1 
            temp(itemp) = etotgatefield/e2
            gatefield_corr => temp(itemp)  
         END IF

         CALL  qexsd_init_total_energy(output%total_energy, etot/e2, eband/e2, ehart/e2, vtxc/e2, &
                                       etxc/e2, ewld/e2, degauss_, demet_, efield_corr, potstat_corr,&
                                       gatefield_corr, DISPERSION_CONTRIBUTION = vdw_term_pt) 
         !
         NULLIFY(degauss_, demet_, efield_corr, potstat_corr, gatefield_corr)
         itemp = 0
          !
!---------------------------------------------------------------------------------------------
! ... FORCES
!----------------------------------------------------------------------------------------------
         !
         IF ( lforce ) THEN 
            output%forces_ispresent = .TRUE.
            CALL qexsd_init_forces(output%forces,nat,force,lforce)
         ELSE 
            output%forces_ispresent = .FALSE.
            output%forces%lwrite = .FALSE.  
         END IF 
         !
!------------------------------------------------------------------------------------------------
! ... STRESS 
!------------------------------------------------------------------------------------------------
         IF ( lstres) THEN
            output%stress_ispresent=.TRUE.
            CALL qexsd_init_stress(output%stress, sigma, lstres ) 
         ELSE 
            output%stress_ispresent=.FALSE.
            output%stress%lwrite=.FALSE.
         END IF
!-------------------------------------------------------------------------------------------------
! ... ELECTRIC FIELD
!-------------------------------------------------------------------------------------------------
         output%electric_field_ispresent = ( gate .OR. lelfield .OR. lberry .OR. tefield ) 

         IF ( gate ) THEN 
            CALL qexsd_init_gate_info(gate_info_temp,"gateInfo", etotgatefield/e2, zgate, nelec, &
                   alat, at, bg, zv, ityp) 
            gate_info_ptr => gate_info_temp    
         END IF             
         IF ( lelfield ) THEN
            ALLOCATE (bp_el_pol(2), bp_ion_pol(3) )
            bp_el_pol = el_pol 
            bp_ion_pol(1:3) = ion_pol(1:3)
         END IF
         IF ( tefield .AND. dipfield) THEN 
            CALL qexsd_init_dipole_info(dipol_obj, el_dipole, ion_dipole, edir, eamp, &
                                  emaxpos, eopreg )  
            dipol_ptr => dipol_obj
         END IF
         IF ( lberry ) bp_obj_ptr => qexsd_bp_obj
         IF (output%electric_field_ispresent) &
            CALL qexsd_init_outputElectricField(output%electric_field, lelfield, tefield, dipfield, &
                 lberry, BP_OBJ = bp_obj_ptr, EL_POL = bp_el_pol, ION_POL = bp_ion_pol,          &
                 GATEINFO = gate_info_ptr, DIPOLE_OBJ =  dipol_ptr) 
         ! 
         IF (ASSOCIATED(gate_info_ptr)) THEN 
            CALL qes_reset (gate_info_ptr)
            NULLIFY(gate_info_ptr)
         ENDIF
         IF (ASSOCIATED (dipol_ptr) ) THEN
            CALL qes_reset (dipol_ptr)
            NULLIFY(dipol_ptr)
         ENDIF
         NULLIFY ( bp_obj_ptr) 
!------------------------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
 10      CONTINUE
         !
         CALL qes_write (qexsd_xf,output)
         CALL qes_reset (output) 
         !
!-------------------------------------------------------------------------------
! ... CLOSING
!-------------------------------------------------------------------------------
         !
         CALL qexsd_closeschema()
         !
      END IF
      DEALLOCATE (ngk_g)
      !
      RETURN
       !
    CONTAINS
       SUBROUTINE check_and_allocate(alloc, mydata)
          IMPLICIT NONE
          REAL(DP),ALLOCATABLE  :: alloc(:) 
          REAL(DP)              :: mydata(:)  
          IF ( ANY(mydata(1:nsp) /= 0.0_DP)) THEN 
             ALLOCATE(alloc(nsp)) 
             alloc(1:nsp) = mydata(1:nsp) 
          END IF 
          RETURN
       END SUBROUTINE check_and_allocate 
    END SUBROUTINE pw_write_schema
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_binaries( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_sum, mp_max
      USE io_base,              ONLY : write_wfc
      USE io_files,             ONLY : iunwfc, nwordwfc
      USE cell_base,            ONLY : tpiba, alat, bg
      USE control_flags,        ONLY : gamma_only, smallmem
      USE gvect,                ONLY : ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol

      USE buffers,              ONLY : get_buffer
      USE wavefunctions, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, igk_k, wk
      USE gvect,                ONLY : ngm, ngm_g, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,                ONLY : ngms_g, dual
      USE wvfct,                ONLY : npwx, et, wg, nbnd
      USE lsda_mod,             ONLY : nspin, isk, lsda
      USE mp_pools,             ONLY : intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : my_bgrp_id, root_bgrp, intra_bgrp_comm,&
                                       root_bgrp_id, nbgrp
      !
      IMPLICIT NONE
      !
      INTEGER               :: i, ig, ngg, ipol, ispin
      INTEGER               :: ik, ik_g, ike, iks, npw_g, npwx_g
      INTEGER, EXTERNAL     :: global_kpoint_index
      INTEGER,  ALLOCATABLE :: ngk_g(:), mill_k(:,:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
      CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=320)    :: filename
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix
      !
      ! ... write wavefunctions and k+G vectors
      !
      iks = global_kpoint_index (nkstot, 1)
      ike = iks + nks - 1
      !
      ! ... ngk_g: global number of k+G vectors for all k points
      !
      ALLOCATE( ngk_g( nkstot ) )
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      CALL mp_sum( ngk_g, inter_pool_comm)
      CALL mp_sum( ngk_g, intra_pool_comm)
      ngk_g = ngk_g / nbgrp
      !
      ! ... npwx_g: maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! ... The igk_l2g array yields the correspondence between the
      ! ... local k+G index and the global G index
      !
      ALLOCATE ( igk_l2g( npwx ) )
      !
      ! ... the igk_l2g_kdip local-to-global map yields the correspondence
      ! ... between the global order of k+G and the local index for k+G.
      !
      ALLOCATE ( igk_l2g_kdip( npwx ) )
      !
      ALLOCATE ( mill_k( 3, npwx ) )
      !
      k_points_loop: DO ik = 1, nks
         !
         ! ik_g is the index of k-point ik in the global list
         !
         ik_g = ik + iks - 1
         !
         ! ... Compute the igk_l2g array from previously computed arrays
         ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
         !
         igk_l2g = 0
         DO ig = 1, ngk (ik)
            igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
         END DO
         !
         ! ... npw_g is the maximum G vector index among all processors
         !
         npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
         CALL mp_max( npw_g, intra_pool_comm )
         !
         igk_l2g_kdip = 0
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik_g), ngk(ik), igk_l2g, &
                              igk_l2g_kdip )
         !
         ! ... mill_k(:,i) contains Miller indices for (k+G)_i
         !
         DO ig = 1, ngk (ik)
            mill_k(:,ig) = mill(:,igk_k(ig,ik))
         END DO
         !
         ! ... read wavefunctions - do not read if already in memory (nsk==1)
         !
         IF ( nks > 1 ) CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
         !
         IF ( nspin == 2 ) THEN
            !
            ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
            !
            ik_g = MOD ( ik_g-1, nkstot/2 ) + 1 
            ispin = isk(ik)
            filename = TRIM(dirname) // 'wfc' // updw(ispin) // &
                 & TRIM(int_to_char(ik_g))
            !
         ELSE
            !
            ispin = 1
            filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_g))
            !
         ENDIF
         !
         ! ... Only the first band group of each pool writes
         ! ... No warranty it works for more than one band group
         !
         IF ( my_bgrp_id == root_bgrp_id ) CALL write_wfc( iunpun, &
              filename, root_bgrp, intra_bgrp_comm, ik_g, tpiba*xk(:,ik), &
              ispin, nspin, evc, npw_g, gamma_only, nbnd, &
              igk_l2g_kdip(:), ngk(ik), tpiba*bg(:,1), tpiba*bg(:,2), &
              tpiba*bg(:,3), mill_k, 1.D0 )
         !
      END DO k_points_loop
      !
      DEALLOCATE ( mill_k )
      DEALLOCATE ( igk_l2g_kdip )
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( ngk_g )
      !
      RETURN
      !
    END SUBROUTINE pw_write_binaries
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gk_l2gmap_kdip( npw_g, ngk_g, ngk, igk_l2g, igk_l2g_kdip, igwk )
      !-----------------------------------------------------------------------
      !
      ! ... This subroutine maps local G+k index to the global G vector index
      ! ... the mapping is used to collect wavefunctions subsets distributed
      ! ... across processors.
      ! ... This map is used to obtained the G+k grids related to each kpt
      !
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE mp,                   ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      ! ... Here the dummy variables
      !
      INTEGER, INTENT(IN)  :: npw_g, ngk_g, ngk
      INTEGER, INTENT(IN)  :: igk_l2g(ngk)
      INTEGER, INTENT(OUT) :: igk_l2g_kdip(ngk)
      INTEGER, OPTIONAL, INTENT(OUT) :: igwk(ngk_g)
      !
      INTEGER, ALLOCATABLE :: igwk_(:), itmp(:), igwk_lup(:)
      INTEGER              :: ig, ig_, ngg
      !
      !
      ALLOCATE( itmp( npw_g ) )
      ALLOCATE( igwk_( ngk_g ) )
      !
      itmp(:)  = 0
      igwk_(:) = 0
      !
      DO ig = 1, ngk
         itmp(igk_l2g(ig)) = igk_l2g(ig)
      END DO
      !
      CALL mp_sum( itmp, intra_bgrp_comm )
      !
      ngg = 0
      DO ig = 1, npw_g
         !
         IF ( itmp(ig) == ig ) THEN
            !
            ngg = ngg + 1
            igwk_(ngg) = ig
            !
         END IF
         !
      END DO
      !
      IF ( ngg /= ngk_g ) &
         CALL errore( 'gk_l2gmap_kdip', 'unexpected dimension in ngg', 1 )
      !
      IF ( PRESENT( igwk ) ) THEN
         !
         igwk(1:ngk_g) = igwk_(1:ngk_g)
         !
      END IF
      !
      ALLOCATE( igwk_lup( npw_g ) )
      !
!$omp parallel private(ig_, ig)
!$omp workshare
      igwk_lup = 0
!$omp end workshare
!$omp do
      DO ig_ = 1, ngk_g
         igwk_lup(igwk_(ig_)) = ig_
      END DO
!$omp end do
!$omp do
      DO ig = 1, ngk
         igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
      END DO
!$omp end do
!$omp end parallel
      !
      DEALLOCATE( igwk_lup )
      !
      DEALLOCATE( itmp, igwk_ )
      !
      RETURN
      !
    END SUBROUTINE gk_l2gmap_kdip

    !------------------------------------------------------------------------
    SUBROUTINE pw_readschema_file(ierr, restart_output, restart_parallel_info, restart_general_info, &
                                  prev_input)
      !------------------------------------------------------------------------
      USE qes_types_module,     ONLY : input_type, output_type, general_info_type, parallel_info_type    
      !
      USE qes_libs_module,      ONLY : qes_write  
      USE FoX_dom,              ONLY : parseFile, item, getElementsByTagname, destroy, nodeList, Node
      USE qes_read_module,      ONLY : qes_read
      IMPLICIT NONE 
      ! 
      INTEGER                                            :: ierr, io_err  
      TYPE( output_type ),OPTIONAL,        INTENT(OUT)   :: restart_output
      TYPE(parallel_info_type),OPTIONAL,   INTENT(OUT)   :: restart_parallel_info
      TYPE(general_info_type ),OPTIONAL,   INTENT(OUT)   :: restart_general_info
      TYPE(input_type),OPTIONAL,           INTENT(OUT)   :: prev_input
      ! 
      TYPE(Node), POINTER     :: root, nodePointer
      TYPE(nodeList),POINTER  :: listPointer
      LOGICAL                 :: found
      CHARACTER(LEN=80)       :: errmsg = ' '
      CHARACTER(LEN=320)      :: filename
      INTEGER,EXTERNAL        :: find_free_unit
      !  
      ! 
      ierr = 0 
      io_err = 0 
      ! 
      iunpun = find_free_unit()
      IF (iunpun .LT. 0 ) &
            CALL errore ("pw_readschema_file", "could not find a free unit to open data-file-schema.xml", 1)
      CALL qexsd_init_schema( iunpun )
      !
      filename = TRIM( tmp_dir ) // TRIM( prefix ) // postfix // TRIM( xmlpun_schema )
      INQUIRE ( file=filename, exist=found )
      IF (.NOT. found ) ierr = ierr + 1
      IF ( ierr /=0 ) THEN
         errmsg='xml data file not found'
         GOTO 100
      END IF
      !
      root => parseFile(filename)
      !
      IF ( PRESENT ( restart_general_info ) ) THEN 
         nodePointer => item ( getElementsByTagname(root, "general_info"),0)
         CALL qes_read( nodePointer, restart_general_info, ierr)
         IF ( ierr /=0 ) THEN
            errmsg='error header of xml data file'
            GOTO 100
         END IF
         ! CALL qes_write_general_info( 82, restart_general_info) 
      END IF 
      ! 
      IF ( PRESENT ( restart_parallel_info ) ) THEN 
         nodePointer => item ( getElementsByTagname(root,"parallel_info"),0)
         CALL qes_read(nodePointer, restart_parallel_info, ierr)
         !
         IF ( ierr /=0) THEN  
            errmsg='error parallel_info  of xsd data file' 
            GOTO 100
         END IF
         ! CALL qes_write_parallel_info ( 82, restart_parallel_info )
      END IF  
      ! 
      IF ( PRESENT ( restart_output ) ) THEN
         nodePointer => item ( getElementsByTagname(root, "output"),0)
         CALL qes_read ( nodePointer, restart_output, ierr ) 
         IF ( ierr /= 0 ) THEN  
            errmsg = 'error output of xsd data file' 
            GOTO 100 
         END IF 
         !
         !CALL qes_write_output ( 82, restart_output ) 
      END IF 
      !
      IF (PRESENT (prev_input)) THEN
         nodePointer => item( getElementsByTagname(root, "input"),0)
         IF ( ASSOCIATED(nodePointer) ) THEN
            CALL qes_read (nodePointer, prev_input, ierr ) 
         ELSE 
            ierr = 5
         END IF
         IF (ierr /= 0 ) THEN
             CALL infomsg ('pw_readschema_file',& 
                            'failed retrieving input info from xml file, please check it')
             IF ( TRIM(prev_input%tagname) == 'input' )  CALL qes_reset (prev_input) 
             ierr = 0 
         END IF
      END IF
      ! 
      CALL destroy(root)       

 100  CALL errore('pw_readschemafile',TRIM(errmsg),ierr)
      !
    END SUBROUTINE pw_readschema_file
    !  
    !------------------------------------------------------------------------
    SUBROUTINE init_vars_from_schema( what, ierr, output_obj, par_info, gen_info, input_obj )
      !------------------------------------------------------------------------
      !
      USE io_rho_xml,           ONLY : read_scf
      USE scf,                  ONLY : rho
      USE lsda_mod,             ONLY : nspin
      USE qes_types_module,     ONLY : input_type, output_type, &
                                       general_info_type, parallel_info_type    
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)           :: what
      TYPE ( output_type), INTENT(IN)        :: output_obj
      TYPE ( parallel_info_type), INTENT(IN) :: par_info
      TYPE ( general_info_type ), INTENT(IN) :: gen_info
      TYPE ( input_type), OPTIONAL, INTENT(IN)         :: input_obj
      INTEGER,INTENT (OUT)                   :: ierr 
      !
      CHARACTER(LEN=256) :: dirname
      LOGICAL            :: lcell, lpw, lions, lspin, linit_mag, &
                            lxc, locc, lbz, lbs, lwfc, lheader,          &
                            lsymm, lrho, lefield, ldim, &
                            lef, lexx, lesm, lpbc, lvalid_input, lalgo, lsymflags 
      !
      LOGICAL            :: found, electric_field_ispresent
      INTEGER            :: tmp
      
      !    
      !
      ierr = 0 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix
      !
      !
      IF ( PRESENT (input_obj) ) THEN 
         lvalid_input = (TRIM(input_obj%tagname) == "input")
      ELSE
         lvalid_input = .FALSE. 
      ENDIF
      !
      !
      ldim    = .FALSE.
      lcell   = .FALSE.
      lpw     = .FALSE.
      lions   = .FALSE.
      lspin   = .FALSE.
      linit_mag = .FALSE.
      lxc     = .FALSE.
      locc    = .FALSE.
      lbz     = .FALSE.
      lbs     = .FALSE.
      lwfc    = .FALSE.
      lsymm   = .FALSE.
      lrho    = .FALSE.
      lefield = .FALSE.
      lef     = .FALSE.
      lexx    = .FALSE.
      lesm    = .FALSE.
      lheader = .FALSE.
      lpbc    = .FALSE.  
      lalgo   = .FALSE. 
      lsymflags =.FALSE. 
      !
     
         
      SELECT CASE( what )
      CASE( 'header' )
         !
         lheader = .TRUE.
         !
      CASE( 'dim' )
         !
         ldim =       .TRUE.
         !
      CASE( 'pseudo' )
         !
         lions = .TRUE.
         !
      CASE( 'config' )
         !
         lcell = .TRUE.
         lions = .TRUE.
         !
      CASE( 'rho' )
         !
         lrho  = .TRUE.
         !
      CASE( 'wave' )
         !
         lpw   = .TRUE.
         lwfc  = .TRUE.
         !
      CASE( 'nowave' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag   = .TRUE.
         lxc     = .TRUE.
         lexx    = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         lalgo   = .TRUE.
         lsymflags = .TRUE.
         !
      CASE( 'all' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag  = .TRUE.
         lxc     = .TRUE.
         lexx    =.TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lwfc    = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         lrho    = .TRUE.
         lpbc    = .TRUE.
         lalgo   = .TRUE. 
         lsymflags = .TRUE. 
         !
      CASE( 'ef' )
         !
         lef        = .TRUE.
         !
      CASE( 'exx' )
         !
         lexx       = .TRUE.
         !
      CASE( 'esm' )
         !
         lesm       = .TRUE.
         !
      CASE( 'boundary_conditions' )  
         !
         lpbc       = .TRUE.
      END SELECT
      !
      !
      IF ( lheader ) THEN 
         CALL readschema_header( gen_info )
      END IF 
      IF ( ldim ) THEN
         ! 
         CALL readschema_dim(par_info, output_obj%atomic_species, &
              output_obj%atomic_structure, output_obj%symmetries, &
              output_obj%basis_set, output_obj%band_structure ) 
         CALL readschema_kdim(output_obj%symmetries, output_obj%band_structure )
         !
      ENDIF
      !
      IF ( lcell ) THEN
         CALL readschema_cell( output_obj%atomic_structure )
      END IF
      !
      IF ( lpw ) THEN
         CALL readschema_planewaves( output_obj%basis_set) 
      END IF
      IF ( lions ) THEN
         CALL readschema_ions( output_obj%atomic_structure, output_obj%atomic_species, dirname)
      END IF
      IF ( lspin ) THEN

         CALL readschema_spin( output_obj%magnetization )
      END IF
      IF (linit_mag) THEN
         CALL readschema_magnetization (  output_obj%band_structure,  output_obj%atomic_species,&
                                          output_obj%magnetization )
      END IF
      IF ( lxc ) THEN
         CALL readschema_xc (  output_obj%atomic_species, output_obj%dft )
      END IF
      IF ( locc ) THEN
         CALL readschema_occupations( output_obj%band_structure )
      END IF
      IF ( lbz ) THEN
         CALL readschema_brillouin_zone( output_obj%symmetries,  output_obj%band_structure )
      END IF
      IF ( lbs ) THEN
         CALL readschema_band_structure( output_obj%band_structure )
      END IF
      IF ( lwfc ) THEN
         !
         IF (output_obj%band_structure%wf_collected)  CALL read_collected_to_evc(dirname ) 
      END IF
      IF ( lsymm ) THEN
         IF ( lvalid_input .and. lsymflags ) THEN 
            CALL readschema_symmetry (  output_obj%symmetries, output_obj%basis_set, input_obj%symmetry_flags )
         ELSE 
            CALL readschema_symmetry( output_obj%symmetries,output_obj%basis_set) 
         ENDIF
      ENDIF
      !
      IF ( lrho ) THEN
         !
         ! ... to read the charge-density we use the routine from io_rho_xml 
         ! ... it also reads ns for ldaU and becsum for PAW
         !
         CALL read_scf( rho, nspin )
         !
      END IF
      IF ( lef ) THEN
               CALL readschema_ef ( output_obj%band_structure) 
         !
      END IF
      ! 
      IF ( lpbc ) THEN
         CALL readschema_outputPBC ( output_obj%boundary_conditions)
      END IF
      !
      IF ( lefield .AND. lvalid_input ) CALL readschema_efield ( input_obj%electric_field ) 
      !
      IF ( lexx .AND. output_obj%dft%hybrid_ispresent  ) CALL readschema_exx ( output_obj%dft%hybrid )
      IF ( lalgo ) CALL readschema_algo(output_obj%algorithmic_info ) 
      !
      RETURN
      !
    END SUBROUTINE init_vars_from_schema
    !-------------------------------------------------------------------------------
    SUBROUTINE readschema_header (gen_info_obj) 
    !-------------------------------------------------------------------------------
    ! 
    USE io_files,            ONLY: qexsd_fmt, qexsd_version, qexsd_init
    USE qes_types_module,    ONLY: general_info_type
    IMPLICIT NONE 
    !
    TYPE (general_info_type ),INTENT(IN)  :: gen_info_obj   
    ! 
    IF ( qexsd_init ) RETURN 
    qexsd_fmt = TRIM (gen_info_obj%xml_format%NAME)
    qexsd_version = TRIM ( gen_info_obj%xml_format%VERSION)
    qexsd_init = .TRUE. 
    !
    END SUBROUTINE readschema_header 
    ! 
    !--------------------------------------------------------------------------
    SUBROUTINE readschema_dim(par_info_obj, atomic_species, atomic_structure, &
         symmetries, basis_set, band_structure ) 
      !
    USE constants,        ONLY : e2
    USE ions_base,        ONLY : nat, nsp
    USE symm_base,        ONLY : nsym
    USE gvect,            ONLY : ngm_g, ecutrho
    USE fft_base,         ONLY : dfftp
    USE gvecs,            ONLY : ngms_g, dual
    USE fft_base,         ONLY : dffts
    USE lsda_mod,         ONLY : lsda
    USE noncollin_module, ONLY : noncolin
    USE klist,            ONLY : nkstot, nelec
    USE wvfct,            ONLY : nbnd, npwx
    USE gvecw,            ONLY : ecutwfc
    USE control_flags,    ONLY : gamma_only
    USE mp_global,        ONLY : nproc_file, nproc_pool_file, &
                                 nproc_image_file, ntask_groups_file, &
                                 nproc_bgrp_file, nproc_ortho_file
    !
    USE qes_types_module, ONLY : parallel_info_type, atomic_species_type, atomic_structure_type, &
                                 symmetries_type, basis_set_type, band_structure_type, input_type  
    IMPLICIT NONE 
    !
    TYPE ( parallel_info_type ),INTENT(IN)     :: par_info_obj
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    TYPE ( symmetries_type ),INTENT(IN)        :: symmetries
    TYPE ( basis_set_type ),INTENT(IN)         :: basis_set
    TYPE ( band_structure_type ),INTENT(IN)    :: band_structure 
    ! 
    INTEGER                                    :: npwx_
    CALL readschema_cell ( atomic_structure ) 
    ! 
    !---------------------------------------------------------------------
    !                                       PARALLEL  DIM 
    !----------------------------------------------------------------------
    nproc_file = par_info_obj%nprocs
    nproc_pool_file = nproc_file/par_info_obj%npool
    nproc_image_file = nproc_file 
    ntask_groups_file = par_info_obj%ntasks
    nproc_bgrp_file = nproc_image_file / par_info_obj%npool / par_info_obj%nbgrp 
    nproc_ortho_file = par_info_obj%ndiag
    !---------------------------------------------------------------------------
    !                                      ATOMS AND SPECIES 
    !--------------------------------------------------------------------------
    nsp = atomic_species%ntyp
    nat = atomic_structure%nat 
    !                                         SIMMETRIES 
    nsym = symmetries%nsym
    !-----------------------------------------------------------------------------
    !                                          BASIS SET 
    !-----------------------------------------------------------------------------
    ecutwfc = basis_set%ecutwfc*e2
    ecutrho = basis_set%ecutrho*e2
    dual = ecutrho/ecutwfc
    npwx_ = basis_set%npwx
    gamma_only= basis_set%gamma_only
    dfftp%nr1 = basis_set%fft_grid%nr1
    dfftp%nr2 = basis_set%fft_grid%nr2          
    dfftp%nr3 = basis_set%fft_grid%nr3
    dffts%nr1 = basis_set%fft_smooth%nr1
    dffts%nr2 = basis_set%fft_smooth%nr2
    dffts%nr3 = basis_set%fft_smooth%nr3
    ngm_g     = basis_set%ngm
    ngms_g    = basis_set%ngms
    !-------------------------------------------------------------------------
    !                                    BAND STRUCTURE  
    !-------------------------------------------------------------------------
    lsda  =    band_structure%lsda
    noncolin = band_structure%noncolin
    nelec =    band_structure%nelec
    nkstot =   band_structure%nks  
    IF (band_structure%nbnd_ispresent) THEN
       nbnd = band_structure%nbnd
    ELSE IF ( band_structure%nbnd_up_ispresent .AND. band_structure%nbnd_dw_ispresent) THEN
       nbnd = ( band_structure%nbnd_up + band_structure%nbnd_dw )
    ELSE 
       CALL errore('init_vars_from_schema: check xml file !!', &
                   'nbnd or nbnd_up+nbnd_dw are missing in band_structure element', 1)
    END IF     
    IF ( lsda ) THEN
       nkstot = nkstot * 2 
       nbnd   = nbnd / 2
    END IF
    END SUBROUTINE readschema_dim
    !
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_cell(atomic_structure )
    !-----------------------------------------------------------------------
    !
    USE constants,         ONLY : pi,tpi
    USE cell_base,         ONLY : ibrav, alat, at, bg, celldm
    USE cell_base,         ONLY : tpiba, tpiba2, omega
    USE qes_types_module,  ONLY : atomic_structure_type
    !
    IMPLICIT NONE 
    ! 
    TYPE ( atomic_structure_type ),INTENT(IN) :: atomic_structure 
    !
    alat = atomic_structure%alat 
    IF ( atomic_structure%bravais_index_ispresent ) THEN 
       ibrav = atomic_structure%bravais_index 
    ELSE 
       ibrav = 0
    END IF
    at(:,1) =  atomic_structure%cell%a1
    at(:,2) =  atomic_structure%cell%a2
    at(:,3) =  atomic_structure%cell%a3
    !! crystal axis are brought into "alat" units
    at = at / alat
    !
    !! if ibrav is present, cell parameters were computed by subroutine
    !! "latgen" using ibrav and celldm parameters: recalculate celldm
    !
    CALL at2celldm (ibrav,alat,at(:,1),at(:,2),at(:,3),celldm)
    !
    tpiba = tpi/alat
    tpiba2= tpiba**2
    CALL volume (alat,at(:,1),at(:,2),at(:,3),omega)
    CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )

    END SUBROUTINE readschema_cell
    ! 
    !------------------------------------------------------------------------
    SUBROUTINE readschema_ions( atomic_structure, atomic_species, dirname ) 
    !------------------------------------------------------------------------
    ! 
    USE ions_base, ONLY : nat, nsp, ityp, amass, atm, tau
    USE cell_base, ONLY : alat
    USE io_files,  ONLY : psfile, pseudo_dir, pseudo_dir_cur
    USE qes_types_module, ONLY: atomic_structure_type, atomic_species_type, input_type 
    ! 
    IMPLICIT NONE 
    ! 
    TYPE ( atomic_structure_type ),INTENT(IN) :: atomic_structure
    TYPE ( atomic_species_type ),INTENT(IN)   :: atomic_species  
    CHARACTER(LEN=*), INTENT(IN)              :: dirname
    ! 
    INTEGER                                   :: iat, isp, idx
    CHARACTER(LEN = 3 ),ALLOCATABLE           :: symbols(:) 
    ! 
    nat = atomic_structure%nat
    nsp = atomic_species%ntyp
    ALLOCATE ( symbols(nat) ) 
    DO isp = 1, nsp 
       amass(isp) = 0.d0 
       IF (atomic_species%species(isp)%mass_ispresent) amass(isp) = atomic_species%species(isp)%mass
       atm(isp) = TRIM ( atomic_species%species(isp)%name )
       psfile(isp) = TRIM ( atomic_species%species(isp)%pseudo_file) 
    END DO 
    !
    loop_on_atoms:DO iat = 1, nat
       idx = atomic_structure%atomic_positions%atom(iat)%index
       tau(:,idx) = atomic_structure%atomic_positions%atom(iat)%atom 
       symbols(idx)  = TRIM ( atomic_structure%atomic_positions%atom(idx)%name ) 
       loop_on_species:DO isp = 1, nsp
          IF ( TRIM(symbols(idx)) == TRIM (atm(isp))) THEN 
             ityp(iat) = isp 
             exit loop_on_species
          END IF 
       END  DO loop_on_species
    END DO loop_on_atoms
    
    DEALLOCATE ( symbols ) 
    IF ( atomic_structure%alat_ispresent ) alat = atomic_structure%alat 
    tau(:,1:nat) = tau(:,1:nat)/alat  
    ! 
    ! ... this is where PP files used in the calculation were stored
    !
    pseudo_dir_cur = TRIM(dirname)
    ! 
    ! ... this is where PP files were originally found (if available)
    !
    IF ( atomic_species%pseudo_dir_ispresent) THEN 
       pseudo_dir = TRIM(atomic_species%pseudo_dir)
    ELSE 
       pseudo_dir = pseudo_dir_cur
    END IF
    ! 
    END SUBROUTINE readschema_ions
    !  
    !------------------------------------------------------------------------
    SUBROUTINE readschema_symmetry ( symms_obj, basis_obj, flags_obj  ) 
    !------------------------------------------------------------------------
      ! 
      USE symm_base,       ONLY : nrot, nsym, invsym, s, ft,ftau, irt, t_rev, &
                                 sname, sr, invs, inverse_s, s_axis_to_cart, &
                                 time_reversal, no_t_rev, nosym
      USE control_flags,   ONLY : noinv  
      USE fft_base,        ONLY : dfftp
      USE qes_types_module,ONLY : symmetries_type, symmetry_type, basis_type
      ! 
      IMPLICIT NONE   
      ! 
      TYPE ( symmetries_type )               :: symms_obj 
      TYPE ( basis_set_type )                :: basis_obj
      TYPE ( symmetry_flags_type),OPTIONAL   :: flags_obj
      INTEGER                                :: isym 
      ! 
      IF ( PRESENT(flags_obj) ) THEN 
         noinv = flags_obj%noinv
         nosym = flags_obj%nosym
         no_t_rev = flags_obj%no_t_rev
         time_reversal = time_reversal .AND. .NOT. noinv 
      ENDIF
      !
      nrot = symms_obj%nrot 
      nsym = symms_obj%nsym
      ! 
      !  
      invsym = .FALSE. 
      DO isym = 1, nrot
        s(:,:,isym) = reshape(symms_obj%symmetry(isym)%rotation%matrix, [3,3]) 
        sname(isym) = TRIM ( symms_obj%symmetry(isym)%info%name )  
        IF ( (TRIM(sname(isym)) == "inversion") .AND. (isym .LE. nsym) ) invsym = .TRUE.
        IF ( symms_obj%symmetry(isym)%fractional_translation_ispresent .AND. (isym .LE. nsym) ) THEN
           ft(1:3,isym)  =  symms_obj%symmetry(isym)%fractional_translation(1:3) 
           ftau(1,isym) = NINT( ft(1,isym)*DBLE( basis_obj%fft_grid%nr1 ) )
           ftau(2,isym) = NINT( ft(2,isym)*DBLE( basis_obj%fft_grid%nr2 ) )
           ftau(3,isym) = NINT( ft(3,isym)*DBLE( basis_obj%fft_grid%nr3 ) )
        END IF
        IF ( symms_obj%symmetry(isym)%info%time_reversal_ispresent ) THEN  
           IF (symms_obj%symmetry(isym)%info%time_reversal) THEN 
              t_rev( isym ) = 1
           ELSE
              t_rev( isym ) = 0 
           END IF 
        END IF
        IF ( symms_obj%symmetry(isym)%equivalent_atoms_ispresent .AND. (isym .LE. nsym) )   &
             irt(isym,:) = symms_obj%symmetry(isym)%equivalent_atoms%equivalent_atoms(:)
      END DO
      CALL inverse_s()
      CALL s_axis_to_cart() 
      !
    END SUBROUTINE readschema_symmetry 
    !
    !---------------------------------------------------------------------------
    SUBROUTINE readschema_efield( efield_obj  ) 
    !---------------------------------------------------------------------------
      !       
      USE extfield, ONLY : tefield, dipfield, edir, emaxpos, eopreg, eamp, gate, zgate, &
                           block, block_1, block_2, block_height, relaxz
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( electric_field_type),OPTIONAL, INTENT(IN)    :: efield_obj
      ! 
      !
      tefield = .FALSE. 
      dipfield = .FALSE. 
      IF ( .NOT. PRESENT( efield_obj) ) RETURN 
      IF (TRIM(efield_obj%electric_potential) == 'sawtooth_potential') THEN 
         tefield = .TRUE. 
         IF ( efield_obj%dipole_correction_ispresent ) THEN 
            dipfield = efield_obj%dipole_correction
         ELSE 
            dipfield = .FALSE. 
         END IF 
         IF ( efield_obj%electric_field_direction_ispresent ) THEN 
            edir = efield_obj%electric_field_direction
         ELSE 
            edir = 3 
         END IF
         IF ( efield_obj%potential_max_position_ispresent ) THEN 
            emaxpos = efield_obj%potential_max_position
         ELSE 
            emaxpos = 5d-1
         END IF 
         IF ( efield_obj%potential_decrease_width_ispresent ) THEN 
            eopreg = efield_obj%potential_decrease_width
         ELSE 
            eopreg = 1.d-1
         END IF 
         IF ( efield_obj%electric_field_amplitude_ispresent ) THEN 
            eamp = efield_obj%electric_field_amplitude
         ELSE 
            eamp = 1.d-3
         END IF
         IF (efield_obj%gate_settings_ispresent) THEN 
            gate = efield_obj%gate_settings%use_gate
            IF (efield_obj%gate_settings%zgate_ispresent) zgate     = efield_obj%gate_settings%zgate
            IF (efield_obj%gate_settings%relaxz_ispresent) relaxz   = efield_obj%gate_settings%relaxz
            IF (efield_obj%gate_settings%block_ispresent) block     = efield_obj%gate_settings%block
            IF (efield_obj%gate_settings%block_1_ispresent) block_1 = efield_obj%gate_settings%block_1
            IF (efield_obj%gate_settings%block_2_ispresent) block_2 = efield_obj%gate_settings%block_2
            IF (efield_obj%gate_settings%block_height_ispresent) &
                                                         block_height = efield_obj%gate_settings%block_height
         END IF 
      END IF 
      !
  END SUBROUTINE readschema_efield  
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_planewaves ( basis_set_obj ) 
    !-----------------------------------------------------------------------
    ! 
    USE constants,       ONLY : e2
    USE gvect,           ONLY : ngm_g, ecutrho
    USE gvecs,           ONLY : ngms_g, dual
    USE gvecw,           ONLY : ecutwfc
    USE fft_base,        ONLY : dfftp
    USE fft_base,        ONLY : dffts
    USE wvfct,           ONLY : npwx
    USE control_flags,   ONLY : gamma_only
    USE qes_types_module,ONLY : basis_set_type
    ! 
    IMPLICIT NONE 
    ! 
    TYPE ( basis_set_type )              :: basis_set_obj  
    !
    ecutwfc = basis_set_obj%ecutwfc*e2
    ecutrho = basis_set_obj%ecutrho*e2
    dual = ecutrho/ecutwfc
    !npwx = basis_set_obj%npwx
    gamma_only= basis_set_obj%gamma_only
    dfftp%nr1 = basis_set_obj%fft_grid%nr1
    dfftp%nr2 = basis_set_obj%fft_grid%nr2          
    dfftp%nr3 = basis_set_obj%fft_grid%nr3
    dffts%nr1 = basis_set_obj%fft_smooth%nr1
    dffts%nr2 = basis_set_obj%fft_smooth%nr2
    dffts%nr3 = basis_set_obj%fft_smooth%nr3
    ngm_g     = basis_set_obj%ngm
    ngms_g    = basis_set_obj%ngms
    !
    END SUBROUTINE readschema_planewaves 
    !--------------------------------------------------------------------------
    SUBROUTINE readschema_spin( magnetization_obj) 
    !--------------------------------------------------------------------------
      ! 
      USE spin_orb,         ONLY : lspinorb, domag
      USE lsda_mod,         ONLY : nspin, lsda
      USE noncollin_module, ONLY : noncolin, npol
      USE qes_types_module, ONLY : magnetization_type
      USE symm_base,        ONLY : time_reversal
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( magnetization_type ),INTENT(IN)         :: magnetization_obj
      ! 
      lspinorb = magnetization_obj%spinorbit 
      domag =   magnetization_obj%do_magnetization 
      lsda  =   magnetization_obj%lsda
      noncolin = magnetization_obj%noncolin  
      IF ( noncolin .AND. domag ) time_reversal = .FALSE.
      IF ( lsda ) THEN  
        nspin = 2
        npol = 1
      ELSE IF (noncolin ) THEN 
        nspin = 4
        npol = 2
      ELSE 
        nspin =1
        npol = 1 
      END IF
      ! 
    END SUBROUTINE readschema_spin 
    !
    !-----------------------------------------------------------------------------------------
    SUBROUTINE readschema_magnetization( band_structure_obj, atomic_specs_obj, magnetization_obj ) 
      !---------------------------------------------------------------------------------------
      ! 
      USE klist,            ONLY : two_fermi_energies, nelup, neldw, tot_magnetization
      USE ener,             ONLY : ef_up, ef_dw
      USE lsda_mod,         ONLY : starting_magnetization
      USE noncollin_module, ONLY : angle1, angle2, i_cons, mcons, bfield, &
                                   lambda
      USE electrons_base,   ONLY : set_nelup_neldw
      USE qes_types_module, ONLY : band_structure_type, atomic_species_type, input_type 
      !
      IMPLICIT NONE 
      ! 
      TYPE ( band_structure_type ),INTENT(IN)    :: band_structure_obj
      TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_specs_obj
      TYPE ( magnetization_type ) ,INTENT(IN)    :: magnetization_obj
      !  
      REAL(DP)                   :: tot_mag_, nelec_, theta, phi, fixed_magnetization(3) 
      INTEGER                    :: nsp_, isp
      !
      bfield = 0.d0
      nelec_ = band_structure_obj%nelec
      two_fermi_energies = band_structure_obj%two_fermi_energies_ispresent
      IF (two_fermi_energies) THEN 
         ef_up = band_structure_obj%two_fermi_energies(1)
         ef_dw = band_structure_obj%two_fermi_energies(2) 
         IF (TRIM(band_structure_obj%occupations_kind%occupations) == 'fixed') THEN
            tot_magnetization = magnetization_obj%total
            CALL set_nelup_neldw(tot_magnetization, nelec_, nelup, neldw) 
         END IF 
      END IF 
      nsp_ = atomic_specs_obj%ntyp
      !
      i_cons = 0
      DO isp = 1, nsp_
         IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent) THEN
             starting_magnetization(isp) = atomic_specs_obj%species(isp)%starting_magnetization      
         END IF                                                                                      
         !                                                                                           
         IF ( band_structure_obj%noncolin ) THEN                                                         
            IF (    atomic_specs_obj%species(isp)%spin_teta_ispresent ) THEN 
               theta = atomic_specs_obj%species(isp)%spin_teta 
               angle1(isp) = theta 
            END IF                                                                  
            IF ( atomic_specs_obj%species(isp)%spin_phi_ispresent ) THEN                  
               phi = atomic_specs_obj%species(isp)%spin_phi
               angle2(isp) = phi
            END IF                                                                     
               !                                                                                     
            IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent .AND. &
                                                                              i_cons == 1 ) THEN 
                !            
                mcons(1,isp) = starting_magnetization(isp) * sin( theta ) * cos( phi )
                mcons(2,isp) = starting_magnetization(isp) * sin( theta ) * sin( phi )
                mcons(3,isp) = starting_magnetization(isp) * cos( theta )
            ELSE IF ( i_cons == 2) THEN  
                mcons(3,isp) = cos(theta) 
            END IF
         ELSE IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent .AND. &
                                                                                  i_cons == 1 ) THEN 
            mcons(1,isp) = starting_magnetization(isp)                                               
         END IF                                                                                      
      END DO   
      !
    END SUBROUTINE readschema_magnetization
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_xc ( atomic_specs, dft_obj ) 
    !-----------------------------------------------------------------------
      ! 
      USE funct,     ONLY : enforce_input_dft
      USE ldaU,      ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax, &
                            Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_alpha, &
                            Hubbard_J0, Hubbard_beta, U_projection
      USE kernel_table,     ONLY : vdw_table_name
      USE control_flags,    ONLY : llondon, lxdm, ts_vdw
      USE london_module,    ONLY : scal6, lon_rcut, in_C6
      USE tsvdw_module,     ONLY : vdw_isolated
      USE qes_types_module, ONLY : atomic_species_type, dft_type
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( atomic_species_type )    :: atomic_specs
      TYPE ( dft_type            )    :: dft_obj
      INTEGER                         :: ihub, nsp_, inlc, isp
      ! 
      CHARACTER(LEN = 256 )           :: label
      CHARACTER(LEN = 20  )           :: dft_name
      CHARACTER(LEN =  3 )            :: symbol
      ! 
      nsp_ = atomic_specs%ntyp 
      dft_name = TRIM(dft_obj%functional) 
      CALL enforce_input_dft ( dft_name, .TRUE. ) 
      lda_plus_u = dft_obj%dftU_ispresent 
      IF ( lda_plus_u ) THEN 
         lda_plus_u_kind = dft_obj%dftU%lda_plus_u_kind
         U_projection = TRIM ( dft_obj%dftU%U_projection_type )
         Hubbard_l =-1 
         IF ( dft_obj%dftU%Hubbard_U_ispresent) THEN 
            loop_on_hubbardU:DO ihub =1, dft_obj%dftU%ndim_Hubbard_U
               symbol = TRIM(dft_obj%dftU%Hubbard_U(ihub)%specie)
               label  = TRIM(dft_obj%dftU%Hubbard_U(ihub)%label ) 
               loop_on_speciesU:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_U(isp) = dft_obj%dftU%Hubbard_U(ihub)%HubbardCommon
                     SELECT CASE ( TRIM (label))
                        CASE ( '1s', '2s', '3s', '4s', '5s', '6s', '7s' ) 
                            Hubbard_l(isp) = 0 
                        CASE ( '2p', '3p', '4p', '5p', '6p' ) 
                            Hubbard_l(isp) = 1 
                        CASE ( '3d', '4d', '5d' ) 
                            Hubbard_l( isp ) = 2 
                        CASE ( '4f', '5f' )  
                            Hubbard_l(isp ) = 3
                        CASE  default 
                            IF (Hubbard_U(isp)/=0) &
                              CALL errore ("pw_readschema:", "unrecognized label for Hubbard "//label, 1 ) 
                     END SELECT   
                     EXIT loop_on_speciesU
                  END IF 
                END DO loop_on_speciesU
            END DO loop_on_hubbardU
         END IF 
         
         IF ( dft_obj%dftU%Hubbard_J0_ispresent ) THEN 
           loop_on_hubbardj0:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J0
               symbol = TRIM(dft_obj%dftU%Hubbard_J0(ihub)%specie)
               loop_on_speciesj0:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN
                     Hubbard_J0(isp) = dft_obj%dftU%Hubbard_J0(ihub)%HubbardCommon
                     EXIT loop_on_speciesj0
                  END IF
               END DO loop_on_speciesj0
            END DO loop_on_hubbardj0
         END IF
         IF ( dft_obj%dftU%Hubbard_alpha_ispresent) THEN 
            loop_on_hubbardAlpha:DO ihub =1, dft_obj%dftU%ndim_Hubbard_alpha
               symbol = TRIM(dft_obj%dftU%Hubbard_alpha(ihub)%specie)
               loop_on_speciesAlpha:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_alpha(isp) = dft_obj%dftU%Hubbard_alpha(ihub)%HubbardCommon
                     EXIT loop_on_speciesAlpha
                  END IF 
                END DO loop_on_speciesAlpha
            END DO loop_on_hubbardAlpha
         END IF 
         IF ( dft_obj%dftU%Hubbard_beta_ispresent) THEN 
            loop_on_hubbardBeta:DO ihub =1, dft_obj%dftU%ndim_Hubbard_beta
               symbol = TRIM(dft_obj%dftU%Hubbard_beta(ihub)%specie)
               loop_on_speciesBeta:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_beta(isp) = dft_obj%dftU%Hubbard_beta(ihub)%HubbardCommon
                     EXIT loop_on_speciesBeta
                  END IF 
                END DO loop_on_speciesBeta
            END DO loop_on_hubbardBeta
         END IF 
         IF ( dft_obj%dftU%Hubbard_J_ispresent) THEN 
            loop_on_hubbardJ:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J
               symbol = TRIM(dft_obj%dftU%Hubbard_J(ihub)%specie)
               loop_on_speciesJ:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_J(:,isp) = dft_obj%dftU%Hubbard_J(ihub)%HubbardJ
                     EXIT loop_on_speciesJ
                  END IF 
                END DO loop_on_speciesJ
            END DO loop_on_hubbardJ
         END IF 
         Hubbard_lmax = MAXVAL( Hubbard_l(1:nsp_) )
      END IF 
      ! 
      IF ( dft_obj%vdW_ispresent ) THEN 
           SELECT CASE( TRIM( dft_obj%vdW%vdw_corr ) )
    !
              CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
                !
                llondon= .TRUE.
                ts_vdw= .FALSE.
                lxdm   = .FALSE.
                !
              CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
                !
                llondon= .FALSE.
                ts_vdw= .TRUE.
                lxdm   = .FALSE.
                !
              CASE( 'XDM', 'xdm' )
                 !
                llondon= .FALSE.
                ts_vdw= .FALSE.
                lxdm   = .TRUE.
                !
              CASE DEFAULT
                !
                llondon= .FALSE.
                ts_vdw = .FALSE.
                lxdm   = .FALSE.
                !
           END SELECT
          ! the following lines set vdw_table_name, if not already set before
          ! (the latter option, added by Yang Jiao, is useful for postprocessing)
          ! NOTA BENE: inlc is not used - this part should be simplified and maybe
          ! moved to somewhere else (e.g. a routine setting the default file name) 
           SELECT CASE ( TRIM (dft_obj%vdW%non_local_term))
             CASE ('vdw1')  
                inlc = 1
             CASE ('vdw2') 
                inlc = 2
             CASE ('vv10' ) 
                inlc = 3 
             CASE ( 'vdW-DF-x') 
                inlc = 4
             CASE ( 'vdW-DF-y')
                inlc = 5
             CASE ( 'vdW-DF-z')
                inlc = 6
             CASE default 
                inlc = 0 
          END SELECT
          IF ( vdw_table_name == ' ' ) THEN 
             IF (inlc == 0 ) THEN 
                vdw_table_name = ''
             ELSE IF ( inlc == 3 ) THEN 
                vdw_table_name = 'rVV10_kernel_table'
             ELSE
                vdw_table_name = 'vdW_kernel_table'
             END IF 
          END IF
          !
          IF (dft_obj%vdW%london_s6_ispresent ) THEN 
             scal6 = dft_obj%vdW%london_s6
          END IF 
          IF ( dft_obj%vdW%london_rcut_ispresent ) THEN 
             lon_rcut = dft_obj%vdW%london_rcut
          END IF
          IF ( dft_obj%vdW%london_c6_ispresent ) THEN
             loop_on_londonC6:DO ihub =1, dft_obj%vdW%ndim_london_c6
                symbol = TRIM(dft_obj%vdW%london_c6(ihub)%specie)
                loop_on_speciesC6:DO isp = 1, nsp_
                   IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                      in_C6(isp) = dft_obj%vdW%london_c6(ihub)%HubbardCommon
                      EXIT loop_on_speciesC6
                   END IF
                END DO loop_on_speciesC6
             END DO loop_on_londonC6
          END IF
          !
          IF (dft_obj%vdW%ts_vdW_isolated_ispresent ) THEN 
             vdW_isolated = dft_obj%vdW%ts_vdW_isolated
          END IF 
      END IF 
      !         
    END SUBROUTINE readschema_xc
    !  
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_kdim( symmetries_obj, band_struct_obj )
    !-----------------------------------------------------------------------------------------------------
       !
       USE lsda_mod,         ONLY : lsda
       USE klist,            ONLY : nkstot
       USE symm_base,        ONLY : nrot 
       USE qes_types_module, ONLY : symmetries_type, band_structure_type
       !
       IMPLICIT NONE
       !
       TYPE ( symmetries_type )    ,INTENT(IN)    :: symmetries_obj 
       TYPE ( band_structure_type ),INTENT(IN)    :: band_struct_obj 
       INTEGER                                    :: nks_
       ! 
       nks_ = band_struct_obj%nks
       nkstot = nks_
       IF ( band_struct_obj%lsda ) nkstot = nkstot * 2  
       !
       nrot = symmetries_obj%nrot
       !
    END SUBROUTINE readschema_kdim    
    !
    ! --------- For 2D cutoff: to read the fact that 2D cutoff was used in scf from new xml----------------
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_outputPBC( boundary_conditions_obj )
    !-----------------------------------------------------------------------------------------------------
       !
       USE Coul_cut_2D,       ONLY : do_cutoff_2D
       !
       IMPLICIT NONE
       !
       TYPE ( outputPBC_type ),INTENT(IN)    :: boundary_conditions_obj 
       ! 
       IF ( TRIM(boundary_conditions_obj%assume_isolated) .EQ. "2D" ) THEN
          do_cutoff_2D=.TRUE.  
       ENDIF
       !
    END SUBROUTINE readschema_outputPBC
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_brillouin_zone( symmetries_obj, band_structure )
    !-----------------------------------------------------------------------------------------------------
       !
       USE lsda_mod, ONLY : lsda, isk
       USE klist,    ONLY : nkstot, xk, wk
       USE start_k,  ONLY : nks_start, xk_start, wk_start, &
                              nk1, nk2, nk3, k1, k2, k3 
       USE symm_base,ONLY : nrot, s, sname
       USE qes_types_module, ONLY : k_points_IBZ_type, occupations_type, symmetries_type, band_structure_type
       !
       IMPLICIT NONE
       !
       TYPE ( symmetries_type ),    INTENT(IN)    :: symmetries_obj 
       TYPE ( band_structure_type ),INTENT(IN)    :: band_structure
       INTEGER                                    :: ik, isym, nks_
       ! 
       nks_ = band_structure%nks
       nkstot = nks_
       IF ( band_structure%lsda ) nkstot = nkstot * 2  
       ! 
       ! 
       DO ik = 1, nks_
          xk(:,ik) = band_structure%ks_energies(ik)%k_point%k_point(:) 
       END DO 
       !!  during lsda computations pw uses, for each k-point in the mesh, a distinct 
       !!  k_point variable for the two spin channels, while in 
       !!  the xml file only one k_point is present
       IF ( band_structure%lsda ) THEN
          DO ik = 1, nks_
             xk(:,nks_+ik) = band_structure%ks_energies(ik)%k_point%k_point(:) 
             isk(ik) = 1
             isk(ik+nks_) = 2
          END DO
       END IF   
       !   
       IF ( band_structure%starting_k_points%monkhorst_pack_ispresent ) THEN 
          nks_start = 0 
          nk1 = band_structure%starting_k_points%monkhorst_pack%nk1 
          nk2 = band_structure%starting_k_points%monkhorst_pack%nk2
          nk3 = band_structure%starting_k_points%monkhorst_pack%nk3 
           k1 = band_structure%starting_k_points%monkhorst_pack%k1
           k2 = band_structure%starting_k_points%monkhorst_pack%k2
           k3 = band_structure%starting_k_points%monkhorst_pack%k3
       ELSE IF (band_structure%starting_k_points%nk_ispresent ) THEN 
           nks_start = band_structure%starting_k_points%nk
           IF ( nks_start > 0 ) THEN 
              IF ( .NOT. ALLOCATED(xk_start) ) ALLOCATE (xk_start(3,nks_start))
              IF ( .NOT. ALLOCATED(wk_start) ) ALLOCATE (wk_start(nks_start))
              IF ( nks_start == size( band_structure%starting_k_points%k_point ) ) THEN 
                 DO ik =1, nks_start
                    xk_start(:,ik) = band_structure%starting_k_points%k_point(ik)%k_point(:) 
                    IF ( band_structure%starting_k_points%k_point(ik)%weight_ispresent) THEN 
                        wk_start(ik) = band_structure%starting_k_points%k_point(ik)%weight 
                    ELSE 
                        wk_start(ik) = 0.d0
                    END IF 
                 END DO
              ELSE
                 CALL infomsg ( "pw_readschema: ", &
                                "actual number of start kpoint not equal to nks_start, set nks_start=0")  
                 nks_start = 0 
              END IF
           END IF
       ELSE 
           CALL errore ("pw_readschema: ", &
                        " no information found for initializing brillouin zone information", 1)
       END IF  
       ! 
       nrot = symmetries_obj%nrot
       DO isym =1, symmetries_obj%ndim_symmetry
          s(:,:,isym)     = reshape(symmetries_obj%symmetry(isym)%rotation%matrix, [3,3])
          sname(isym) = TRIM ( symmetries_obj%symmetry(isym)%info%name) 
       END DO 
       !
    END SUBROUTINE readschema_brillouin_zone     
    !--------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_occupations( band_struct_obj ) 
      !------------------------------------------------------------------------------------------------
      ! 
      USE lsda_mod,         ONLY : lsda, nspin
      USE fixed_occ,        ONLY : tfixed_occ, f_inp
      USE ktetra,           ONLY : ntetra, tetra_type
      USE klist,            ONLY : ltetra, lgauss, ngauss, degauss, smearing
      USE electrons_base,   ONLY : nupdwn 
      USE wvfct,            ONLY : nbnd
      USE input_parameters, ONLY : input_parameters_occupations => occupations
      USE qes_types_module, ONLY : input_type, band_structure_type
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( band_structure_type ),INTENT(IN)     :: band_struct_obj 
      INTEGER                                     :: ispin, nk1, nk2, nk3, aux_dim1, aux_dim2 
      ! 
      lsda= band_struct_obj%lsda
      nbnd = band_struct_obj%nbnd
      IF ( band_struct_obj%nbnd_up_ispresent ) nupdwn(1) = band_struct_obj%nbnd_up
      IF ( band_struct_obj%nbnd_dw_ispresent ) nupdwn(2) = band_struct_obj%nbnd_dw 
      IF ( lsda )  THEN 
         nspin = 2
         nbnd = nbnd / 2
      ELSE IF ( band_struct_obj%noncolin) THEN 
         nspin = 4 
      ELSE 
         nspin = 1 
      END IF 
      !
      lgauss = .FALSE. 
      ltetra = .FALSE. 
      tetra_type = 0
      ngauss = 0
      input_parameters_occupations = TRIM ( band_struct_obj%occupations_kind%occupations ) 
      IF (TRIM(input_parameters_occupations) == 'tetrahedra' ) THEN 
        ltetra = .TRUE. 
        nk1 = band_struct_obj%starting_k_points%monkhorst_pack%nk1
        nk2 = band_struct_obj%starting_k_points%monkhorst_pack%nk2
        nk3 = band_struct_obj%starting_k_points%monkhorst_pack%nk3
        ntetra = 6* nk1 * nk2 * nk3 
      ELSE IF (TRIM(input_parameters_occupations) == 'tetrahedra_lin' .OR. &
               TRIM(input_parameters_occupations) == 'tetrahedra-lin' ) THEN
        ltetra = .TRUE. 
        nk1 = band_struct_obj%starting_k_points%monkhorst_pack%nk1
        nk2 = band_struct_obj%starting_k_points%monkhorst_pack%nk2
        nk3 = band_struct_obj%starting_k_points%monkhorst_pack%nk3
        tetra_type = 1
        ntetra = 6* nk1 * nk2 * nk3 
      ELSE IF (TRIM(input_parameters_occupations) == 'tetrahedra_opt' .OR. &
               TRIM(input_parameters_occupations) == 'tetrahedra-opt' ) THEN 
        ltetra = .TRUE. 
        nk1 = band_struct_obj%starting_k_points%monkhorst_pack%nk1
        nk2 = band_struct_obj%starting_k_points%monkhorst_pack%nk2
        nk3 = band_struct_obj%starting_k_points%monkhorst_pack%nk3
        tetra_type = 2
        ntetra = 6* nk1 * nk2 * nk3 
      ELSE IF ( TRIM (input_parameters_occupations) == 'smearing') THEN 
        lgauss = .TRUE.  
        degauss = band_struct_obj%smearing%degauss
        SELECT CASE ( TRIM( band_struct_obj%smearing%smearing ) )
           CASE ( 'gaussian', 'gauss', 'Gaussian', 'Gauss' )
             ngauss = 0
             smearing  = 'gaussian'
           CASE ( 'methfessel-paxton', 'm-p', 'mp', 'Methfessel-Paxton', 'M-P', 'MP' )
             ngauss = 1
             smearing = 'mp'
           CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt', 'M-V', 'MV')
             ngauss = -1
             smearing  = 'mv'
           CASE ( 'fermi-dirac', 'f-d', 'fd', 'Fermi-Dirac', 'F-D', 'FD')
             ngauss = -99
             smearing = 'fd'
        END SELECT
      END IF       
     !
    END SUBROUTINE readschema_occupations
 !
    !------------------------------------------------------------------------
    SUBROUTINE readschema_band_structure( band_struct_obj )
      !------------------------------------------------------------------------
      !
      USE constants,     ONLY : e2
      USE basis,    ONLY : natomwfc
      USE lsda_mod, ONLY : lsda, isk
      USE klist,    ONLY : nkstot, wk, nelec
      USE wvfct,    ONLY : et, wg, nbnd
      USE ener,     ONLY : ef, ef_up, ef_dw
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type)         :: band_struct_obj
      INTEGER                             :: ik, nbnd_, nbnd_up_, nbnd_dw_
      ! 
      !! left here to write bw compatible xml
      lsda = band_struct_obj%lsda
      nkstot = band_struct_obj%nks 
      IF ( lsda) THEN 
         IF (band_struct_obj%nbnd_ispresent) THEN 
            nbnd  = band_struct_obj%nbnd / 2
         ELSE IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent ) THEN 
            nbnd = (band_struct_obj%nbnd_up + band_struct_obj%nbnd_dw)/2 
         ELSE 
            CALL errore ('init_vars_from_schema: ','band_structure xml element nbnd and nbnd_up+nbnd_dw missing', 1)  
         END IF 
         nkstot = nkstot * 2 
         isk(1:nkstot/2) = 1
         isk(nkstot/2+1:nkstot) = 2 
      ELSE 
         isk(1:nkstot)   = 1 
      END IF 
      ! 
      nelec = band_struct_obj%nelec
      natomwfc = band_struct_obj%num_of_atomic_wfc
      IF ( band_struct_obj%fermi_energy_ispresent) THEN 
         ef = band_struct_obj%fermi_energy*e2 
      ELSE IF ( band_struct_obj%two_fermi_energies_ispresent ) THEN 
         ef = 0.d0 
         ef_up = band_struct_obj%two_fermi_energies(1)*e2
         ef_dw = band_struct_obj%two_fermi_energies(2)*e2
      ELSE 
         ef = 0.d0
         ef_up = 0.d0
         ef_dw = 0.d0
      END IF 
      DO ik =1, band_struct_obj%ndim_ks_energies
         IF ( band_struct_obj%lsda) THEN
            IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent) THEN
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%nbnd_dw 
            ELSE IF ( band_struct_obj%nbnd_up_ispresent ) THEN 
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%eigenvalues%size - nbnd_up_
            ELSE IF ( band_struct_obj%nbnd_dw_ispresent ) THEN 
               nbnd_dw_ = band_struct_obj%nbnd_dw
               nbnd_up_ = band_struct_obj%ks_energies(ik)%eigenvalues%size - nbnd_dw_ 
            ELSE 
               nbnd_up_ = band_struct_obj%ks_energies(ik)%eigenvalues%size/2  
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%eigenvalues%size/2
            END IF
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            wk( ik + band_struct_obj%ndim_ks_energies ) = wk(ik) 
            et(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues%vector(1:nbnd_up_)*e2
            et(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                             band_struct_obj%ks_energies(ik)%eigenvalues%vector(nbnd_up_+1:nbnd_up_+nbnd_dw_)*e2
            wg(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%occupations%vector(1:nbnd_up_)*wk(ik)
            wg(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                             band_struct_obj%ks_energies(ik)%occupations%vector(nbnd_up_+1:nbnd_up_+nbnd_dw_)*wk(ik)
         ELSE 
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            nbnd_ = band_struct_obj%ks_energies(ik)%eigenvalues%size
            et (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues%vector(1:nbnd_)*e2
            wg (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%occupations%vector(1:nbnd_)*wk(ik)
         END IF  
      END DO 
    END SUBROUTINE readschema_band_structure 
    !
    !--------------------------------------------------------------------------
    SUBROUTINE readschema_algo(algo_obj) 
       USE control_flags, ONLY: tqr 
       USE realus,        ONLY: real_space 
       IMPLICIT NONE 
       TYPE(algorithmic_info_type),INTENT(IN)   :: algo_obj
       tqr = algo_obj%real_space_q 
       real_space = algo_obj%real_space_beta 
    END SUBROUTINE readschema_algo 

    ! 
    !------------------------------------------------------------------------
    SUBROUTINE read_collected_to_evc( dirname )
      !------------------------------------------------------------------------
      !
      ! ... This routines reads wavefunctions from the new file format and
      ! ... writes them into the old format
      !
      USE control_flags,        ONLY : gamma_only
      USE lsda_mod,             ONLY : nspin, isk
      USE klist,                ONLY : nkstot, wk, nks, xk, ngk, igk_k
      USE wvfct,                ONLY : npwx, g2kin, et, wg, nbnd
      USE wavefunctions, ONLY : evc
      USE io_files,             ONLY : nwordwfc, iunwfc
      USE buffers,              ONLY : save_buffer
      USE gvect,                ONLY : ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol
      USE mp_bands,             ONLY : nbgrp, root_bgrp, intra_bgrp_comm
      USE mp_pools,             ONLY : me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp,                   ONLY : mp_sum, mp_max
      USE io_base,              ONLY : read_wfc
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      !
      CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
      CHARACTER(LEN=320)   :: filename, msg
      INTEGER              :: i, ik, ik_g, ig, ipol, ik_s
      INTEGER              :: npol_, npwx_g, nbnd_
      INTEGER              :: nupdwn(2), ike, iks, npw_g, ispin
      INTEGER, EXTERNAL    :: global_kpoint_index
      INTEGER, ALLOCATABLE :: ngk_g(:), mill_k(:,:)
      INTEGER, ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
      LOGICAL              :: opnd, ionode_k
      REAL(DP)             :: scalef, xk_(3), b1(3), b2(3), b3(3)
 
      !
      iks = global_kpoint_index (nkstot, 1)
      ike = iks + nks - 1
      !
      ! ... ngk_g: global number of k+G vectors for all k points
      !
      ALLOCATE( ngk_g( nkstot ) )
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      CALL mp_sum( ngk_g, inter_pool_comm)
      CALL mp_sum( ngk_g, intra_pool_comm)
      ngk_g = ngk_g / nbgrp
      !
      ! ... npwx_g: maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! ... the root processor of each pool reads
      !
      ionode_k = (me_pool == root_pool)
      !
      ! ... The igk_l2g array yields the correspondence between the
      ! ... local k+G index and the global G index
      !
      ALLOCATE ( igk_l2g( npwx ) )
      !
      ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
      !
      ALLOCATE ( igk_l2g_kdip( npwx ) )
      !
      ALLOCATE( mill_k ( 3,npwx ) )
      !
      k_points_loop: DO ik = 1, nks
         !
         ! index of k-point ik in the global list
         !
         ik_g = ik + iks - 1
         !
         ! ... Compute the igk_l2g array from previously computed arrays
         ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
         !
         igk_l2g = 0
         DO ig = 1, ngk(ik)
            igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
         END DO
         !
         ! ... npw_g: the maximum G vector index among all processors
         !
         npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
         CALL mp_max( npw_g, intra_pool_comm )
         !
         igk_l2g_kdip = 0
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik_g), ngk(ik), igk_l2g, &
                              igk_l2g_kdip )
         !
         evc=(0.0_DP, 0.0_DP)
         !
         IF ( nspin == 2 ) THEN
            !
            ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
            !
            ik_g = MOD ( ik_g-1, nkstot/2 ) + 1 
            ispin = isk(ik)
            filename = TRIM(dirname) // 'wfc' // updw(ispin) // &
                 & TRIM(int_to_char(ik_g))
            !
         ELSE
            !
            filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_g))
            !
         ENDIF
         !
         CALL read_wfc( iunpun, filename, root_bgrp, intra_bgrp_comm, &
              ik_g, xk_, ispin, npol_, evc, npw_g, gamma_only, nbnd_, &
              igk_l2g_kdip(:), ngk(ik), b1, b2, b3, mill_k, scalef )
         !
         ! ... here one should check for consistency between what is read
         ! ... and what is expected
         !
         IF ( nbnd_ < nbnd ) THEN
            WRITE (msg,'("The number of bands for this run is",I6,", but only",&
                 & I6," bands were read from file")')  nbnd, nbnd_  
            CALL errore ('pw_restart - read_collected_to_evc', msg, 1 )
         END IF
         CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
         ! 
      END DO k_points_loop
      !
      DEALLOCATE ( mill_k )
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( igk_l2g_kdip )
      !
      RETURN
      !
    END SUBROUTINE read_collected_to_evc
    !
    !----------------------------------------------------------------------------------------
    SUBROUTINE readschema_ef ( band_struct_obj )
    !----------------------------------------------------------------------------------------
       !
       USE constants, ONLY        : e2
       USE ener,  ONLY            : ef, ef_up, ef_dw
       USE klist, ONLY            : two_fermi_energies, nelec
       USE qes_types_module, ONLY : band_structure_type 
       ! 
       IMPLICIT NONE 
       ! 
       TYPE ( band_structure_type ),INTENT(IN)      :: band_struct_obj 
       ! 
       two_fermi_energies = band_struct_obj%two_fermi_energies_ispresent 
       nelec = band_struct_obj%nelec
       IF ( two_fermi_energies) THEN 
          ef_up = band_struct_obj%two_fermi_energies(1)*e2
          ef_dw = band_struct_obj%two_fermi_energies(2)*e2
       ELSE IF ( band_struct_obj%fermi_energy_ispresent ) THEN 
          ef = band_struct_obj%fermi_energy*e2
       END IF 
    END SUBROUTINE readschema_ef 
    !------------------------------------------------------------------------
    SUBROUTINE readschema_exx ( hybrid_obj) 
    !------------------------------------------------------------------------
      ! 
      USE constants,            ONLY : e2
      USE funct,                ONLY : set_exx_fraction, set_screening_parameter, &
                                      set_gau_parameter, enforce_input_dft, start_exx
      USE exx_base,             ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut
      USE exx,                  ONLY : ecutfock, local_thr
      ! 
      USE  qes_types_module,   ONLY : hybrid_type 
      IMPLICIT NONE
      ! 
      TYPE ( hybrid_type), INTENT(IN)               :: hybrid_obj 
      ! 
      x_gamma_extrapolation = hybrid_obj%x_gamma_extrapolation
      nq1 = hybrid_obj%qpoint_grid%nqx1
      nq2 = hybrid_obj%qpoint_grid%nqx2
      nq3 = hybrid_obj%qpoint_grid%nqx3
      CALL set_exx_fraction( hybrid_obj%exx_fraction) 
      CALL set_screening_parameter ( hybrid_obj%screening_parameter)
      exxdiv_treatment = hybrid_obj%exxdiv_treatment 
      ecutvcut = hybrid_obj%ecutvcut*e2
      ecutfock = hybrid_obj%ecutfock*e2
      IF (hybrid_obj%localization_threshold_ispresent) THEN
         local_thr = hybrid_obj%localization_threshold  
      ELSE 
         local_thr = 0._DP 
      END IF 
      CALL start_exx() 
    END SUBROUTINE  readschema_exx 
    !-----------------------------------------------------------------------------------  

  END MODULE pw_restart_new
