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
  ! ... New PWscf I/O using xml schema and (optionally) hdf5 binaries
  ! ... Parallel execution: the xml file is written by one processor only
  ! ... ("ionode_id"), read by all processors ;
  ! ... the wavefunction files are written / read by one processor per pool,
  ! ... collected on / distributed to all other processors in pool
  !
  USE kinds, ONLY: dp
  USE qes_types_module
  USE qes_write_module, ONLY: qes_write
  USE qes_reset_module, ONLY: qes_reset 
  USE qexsd_module, ONLY: qexsd_openschema, qexsd_closeschema, qexsd_xf
  USE qexsd_input,  ONLY: qexsd_input_obj
  USE qexsd_init,   ONLY: qexsd_init_convergence_info, qexsd_init_algorithmic_info,    & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure,      &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure,          &
                          qexsd_init_dipole_info, qexsd_init_total_energy,             &
                          qexsd_init_forces, qexsd_init_stress,                        &
                          qexsd_init_outputElectricField, qexsd_occ_obj,               &
                          qexsd_init_outputPBC, qexsd_init_gate_info, qexsd_init_hybrid,&
                          qexsd_init_dftU, qexsd_init_vdw
  USE io_global, ONLY : ionode, ionode_id
  USE io_files,  ONLY : iunpun, xmlfile
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  PRIVATE
  PUBLIC :: pw_write_schema, pw_write_binaries, read_collected_to_evc
  !
  CONTAINS
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_schema( only_init, wf_collect )
      !------------------------------------------------------------------------
      !
      ! only_init  = T  write only variables that are known after the 
      !                 initial steps of initialization (e.g. structure)
      !            = F  write the complete xml file
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
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, zv, amass
      USE noncollin_module,     ONLY : noncolin, npol
      USE io_files,             ONLY : psfile, pseudo_dir
      USE klist,                ONLY : nks, nkstot, xk, ngk, wk, &
                                       lgauss, ngauss, smearing, degauss, nelec, &
                                       two_fermi_energies, nelup, neldw, tot_charge, ltetra 
      USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE gvect,                ONLY : ngm, ngm_g, g
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
      USE qexsd_init,           ONLY : qexsd_bp_obj, qexsd_start_k_obj
      USE qexsd_input,          ONLY : qexsd_init_k_points_ibz, &
              qexsd_init_occupations, qexsd_init_smearing
      USE fcp_variables,        ONLY : lfcpopt, lfcpdyn, fcp_mu  
      USE control_flags,        ONLY : conv_elec, conv_ions, ldftd3, do_makov_payne 
      USE Coul_cut_2D,          ONLY : do_cutoff_2D 
      USE esm,                  ONLY : do_comp_esm 
      USE martyna_tuckerman,    ONLY : do_comp_mt 
      USE run_info,             ONLY : title
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(IN) :: only_init, wf_collect
      !
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=8)      :: smearing_loc
      CHARACTER(LEN=8), EXTERNAL :: schema_smearing
      INTEGER               :: i, ig, ngg, ipol
      INTEGER               :: npwx_g, ispin, inlc
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      LOGICAL               :: occupations_are_fixed
      INTEGER                  :: iclass, isym, ielem
      CHARACTER(LEN=15)        :: symop_2_class(48)
      LOGICAL                  :: opt_conv_ispresent, dft_is_vdw, empirical_vdw
      INTEGER                  :: n_opt_steps, n_scf_steps_, h_band
      REAL(DP),TARGET                 :: h_energy
      TYPE(gateInfo_type),TARGET      :: gate_info_temp
      TYPE(gateInfo_type),POINTER     :: gate_info_ptr
      TYPE(dipoleOutput_type),TARGET  :: dipol_obj 
      TYPE(dipoleOutput_type),POINTER :: dipol_ptr 
      TYPE(BerryPhaseOutput_type),  POINTER :: bp_obj_ptr
      TYPE(hybrid_type), POINTER            :: hybrid_obj 
      TYPE(vdW_type), POINTER               :: vdw_obj
      TYPE(dftU_type), POINTER              :: dftU_obj 
      REAL(DP), TARGET                      :: lumo_tmp, ef_targ, dispersion_energy_term 
      REAL(DP), POINTER                     :: lumo_energy, ef_point 
      REAL(DP), ALLOCATABLE                 :: ef_updw(:)
      !
      !
      !
      TYPE(output_type) :: output
      REAL(DP),POINTER    :: degauss_, demet_, efield_corr, potstat_corr,  gatefield_corr  
      LOGICAL, POINTER    :: optimization_has_converged 
      LOGICAL, TARGET     :: conv_opt  
      LOGICAL             :: scf_has_converged 
      INTEGER             :: itemp = 1
      REAL(DP),ALLOCATABLE :: london_c6_(:), bp_el_pol(:), bp_ion_pol(:), U_opt(:), J0_opt(:), alpha_opt(:), &
                              J_opt(:,:), beta_opt(:) 
      CHARACTER(LEN=3),ALLOCATABLE :: species_(:)
      CHARACTER(LEN=20),TARGET   :: dft_nonlocc_
      INTEGER,TARGET             :: dftd3_version_
      CHARACTER(LEN=20),TARGET   :: vdw_corr_, pbc_label 
      CHARACTER(LEN=20),POINTER  :: non_local_term_pt, vdw_corr_pt 
      REAL(DP),TARGET            :: temp(20), lond_rcut_, lond_s6_, ts_vdw_econv_thr_, xdm_a1_, xdm_a2_, ectuvcut_,&
                                    scr_par_, loc_thr_  
      REAL(DP),POINTER           :: vdw_term_pt, ts_thr_pt, london_s6_pt, london_rcut_pt, xdm_a1_pt, xdm_a2_pt, &
                                    ts_vdw_econv_thr_pt, ectuvcut_opt, scr_par_opt, loc_thr_p, h_energy_ptr
      LOGICAL,TARGET             :: dftd3_threebody_, ts_vdw_isolated_
      LOGICAL,POINTER            :: ts_isol_pt, dftd3_threebody_pt, ts_vdw_isolated_pt 
      INTEGER,POINTER            :: dftd3_version_pt
      TYPE(smearing_type),TARGET :: smear_obj 
      TYPE(smearing_type),POINTER:: smear_obj_ptr 

      NULLIFY( degauss_, demet_, efield_corr, potstat_corr, gatefield_corr) 
      NULLIFY( gate_info_ptr, dipol_ptr, bp_obj_ptr, hybrid_obj, vdw_obj, dftU_obj, lumo_energy, ef_point)  
      NULLIFY ( optimization_has_converged, non_local_term_pt, vdw_corr_pt, vdw_term_pt, ts_thr_pt, london_s6_pt,  &
                xdm_a1_pt, xdm_a2_pt, ts_vdw_econv_thr_pt, ts_isol_pt, dftd3_threebody_pt, ts_vdw_isolated_pt,     & 
                dftd3_version_pt )
      NULLIFY ( ectuvcut_opt, scr_par_opt, loc_thr_p, h_energy_ptr, smear_obj_ptr) 

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
      ! XML descriptor
      ! 
      IF ( ionode ) THEN  
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
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
                scf_has_converged = conv_elec 
                n_scf_steps_ = n_scf_steps
         END SELECT
         ! 
            call qexsd_init_convergence_info(output%convergence_info,   &
                        SCf_HAS_CONVERGED = scf_has_converged, &
                        OPTIMIZATION_HAS_CONVERGED = optimization_has_converged,& 
                        N_SCF_STEPS = n_scf_steps_, SCF_ERROR=scf_error/e2,&
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
                    ts_vdw_isolated_ = vdw_isolated
                    ts_vdw_isolated_pt => ts_vdw_isolated_
                    ts_vdw_econv_thr_ = vdw_econv_thr
                    ts_vdw_econv_thr_pt => ts_vdw_econv_thr_
                END IF
            ELSE
                vdw_corr_ = 'none'
                vdw_corr_pt => vdw_corr_
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
         IF (ANY([do_makov_payne, do_comp_mt, do_comp_esm, do_cutoff_2D]))  THEN
            output%boundary_conditions_ispresent=.TRUE.
            IF (do_makov_payne) THEN 
               pbc_label = 'makov_payne' 
            ELSE IF ( do_comp_mt) THEN 
               pbc_label = 'martyna_tuckerman' 
            ELSE IF ( do_comp_esm) THEN 
               pbc_label = 'esm' 
            ELSE IF ( do_cutoff_2D) THEN 
               pbc_label = '2D'
            ELSE 
               CALL errore ('pw_restart_new.f90: ', 'internal error line 470', 1) 
            END IF 
            CALL qexsd_init_outputPBC(output%boundary_conditions, TRIM(pbc_label) )  
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
         IF ( only_init ) GO TO 10
         !
         IF ( .NOT. ( lgauss .OR. ltetra )) THEN 
            occupations_are_fixed = .TRUE.
            CALL get_homo_lumo( h_energy, lumo_tmp)
            h_energy = h_energy/e2
            h_energy_ptr => h_energy 
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
         ELSE IF (ltetra .OR. lgauss) THEN  
                ef_targ = ef/e2
                ef_point => ef_targ
         END IF


         IF ( lgauss ) THEN
            IF (TRIM(qexsd_input_obj%tagname) == 'input') THEN 
               smear_obj = qexsd_input_obj%bands%smearing
            ELSE
               smearing_loc = schema_smearing( smearing )
               CALL qexsd_init_smearing(smear_obj, smearing_loc, degauss)
            END IF  
            smear_obj_ptr => smear_obj  
         END IF 
         !  
            
         CALL qexsd_init_band_structure(  output%band_structure,lsda,noncolin,lspinorb, nelec, natomwfc, &
                                 et, wg, nkstot, xk, ngk_g, wk, SMEARING = smear_obj_ptr,  &
                                 STARTING_KPOINTS = qexsd_start_k_obj, OCCUPATIONS_KIND = qexsd_occ_obj, &
                                 WF_COLLECTED = wf_collect, NBND = nbnd, FERMI_ENERGY = ef_point, EF_UPDW = ef_updw,& 
                                 HOMO = h_energy_ptr, LUMO = lumo_energy )
         ! 
         IF (lgauss)  CALL qes_reset (smear_obj)
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
!-------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
 10      CONTINUE
         !
         CALL qexsd_openschema( xmlfile(), iunpun, 'PWSCF', title )
         CALL qes_write (qexsd_xf,output)
         CALL qes_reset (output) 
         CALL qexsd_closeschema()
         !
!-------------------------------------------------------------------------------
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
       !
    END SUBROUTINE pw_write_schema
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_binaries( )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_sum, mp_max
      USE io_base,              ONLY : write_wfc
      USE io_files,             ONLY : restart_dir, iunwfc, nwordwfc
      USE cell_base,            ONLY : tpiba, alat, bg
      USE control_flags,        ONLY : gamma_only, smallmem
      USE gvect,                ONLY : ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol

      USE buffers,              ONLY : get_buffer
      USE wavefunctions, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, igk_k, wk
      USE gvect,                ONLY : ngm, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE wvfct,                ONLY : npwx, et, wg, nbnd
      USE lsda_mod,             ONLY : nspin, isk, lsda
      USE mp_pools,             ONLY : intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : my_bgrp_id, root_bgrp, intra_bgrp_comm,&
                                       root_bgrp_id
      !
      IMPLICIT NONE
      !
      INTEGER               :: i, ig, ngg, ipol, ispin
      INTEGER               :: ik, ik_g, ike, iks, npw_g
      INTEGER, EXTERNAL     :: global_kpoint_index
      INTEGER,  ALLOCATABLE :: ngk_g(:), mill_k(:,:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
      CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=320)    :: filename
      !
      dirname = restart_dir () 
      !
      ! ... write wavefunctions and k+G vectors
      !
      iks = global_kpoint_index (nkstot, 1)
      ike = iks + nks - 1
      !
      ! ... ngk_g: global number of k+G vectors
      !
      ALLOCATE( ngk_g( nks ) )
      ngk_g(1:nks) = ngk(1:nks)
      CALL mp_sum( ngk_g, intra_bgrp_comm)
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
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik), igk_l2g, &
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
      USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm
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
      INTEGER              :: npol_, nbnd_
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
      ALLOCATE( ngk_g( nks ) )
      ngk_g(1:nks) = ngk(1:nks)
      CALL mp_sum( ngk_g, intra_bgrp_comm)
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
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik), igk_l2g, &
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
      DEALLOCATE ( ngk_g )
      !
      RETURN
      !
    END SUBROUTINE read_collected_to_evc
    !
    !------------------------------------------------------------------------
  END MODULE pw_restart_new
