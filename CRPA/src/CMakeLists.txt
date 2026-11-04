set(src_crpa
	crpa_allocate_q.f90
	crpa_bcast_input.f90
	crpa_clean_q.f90
	crpa_close_q.f90
	crpacom.f90
	crpa_compute_perturbation_q.f90
	crpa_copy_original_xk.f90
	crpa_dealloc_1.f90
	crpa_dealloc_2.f90
	crpa_dealloc_q.f90
	crpa_dvpsi_pert.f90
	crpa_generate_grids.f90
	crpa_init.f90
	crpa_init_q.f90
	crpa_load_q.f90
	crpa_openfil_q.f90
	crpa_postprocess_coulomb.f90
	crpa_prepare_perturbations_q.f90
	crpa_prepare_q.f90
	crpa_print_clock.f90
	crpa_psym_dmag.f90
	crpa_psymdvscf.f90
	crpa_q_points.f90
	crpa_readin.f90
	crpa_R_points.f90
	crpa_run_nscf.f90
	crpa_setup_q.f90      
	crpa_setup_wannier.f90
	crpa_skip_nscf_by_rotate.f90
	crpa_solve_linear_system.f90
	crpa_stop_smoothly.f90
	crpa_summary.f90
	crpa_summary_q.f90
	crpa_sym_dmag.f90
	crpa_symdvscf.f90
	crpa_symmetry.f90
	crpa_update_coulomb.f90
	lr_compute_matrix_elements.f90
	crpa_q0correction.f90
	crpa_phescf_w.f90
	crpa_solve_e_w.f90
	crpa_dielec.f90
	crpa_plasma_frequency.f90
	crpa_diel_lanczos.f90
	crpa_allocate_diel_lanczos.f90
	crpa_do_lanczos_epsilon.f90
	crpa_psym_diel.f90
	crpa_calc_chi_epsilon.f90
	crpa_extrapolate_diel.f90
	crpa_restart_diel_lanczos.f90
	crpa_build_freq_grid.f90
	crpa_coul_lanczos.f90
	crpa_allocate_coul_lanczos.f90
	crpa_do_lanczos_coulomb.f90
	crpa_calc_coulomb_lanczos.f90
	crpa_extrapolate_coul.f90
	crpa_diel_lanczos_chain_io.f90
	crpa_coul_lanczos_chain_io.f90
	crpa_restart_coul_lanczos.f90
	crpa_restart_coul_cg.f90
	crpa_average_w.f90
	crpa_run_spectrum_only.f90
	)

qe_add_library(qe_crpa ${src_crpa})
target_link_libraries( qe_crpa
    PRIVATE qe_pw
	    qe_crpa_w90
            qe_lr_modules
            qe_modules
            qe_fftx
            qe_lax
            qe_upflib
            qe_xclib
            qe_phonon_ph)

###########################################################
# crpa.x
###########################################################
set(src_crpa_x crpa_main.f90)
qe_add_executable(qe_crpa_exe ${src_crpa_x})
set_target_properties(qe_crpa_exe PROPERTIES OUTPUT_NAME crpa.x)
target_link_libraries(qe_crpa_exe PRIVATE qe_pw qe_crpa_w90  qe_lax qe_lr_modules qe_modules qe_crpa qe_xclib qe_fftx qe_phonon_ph)
