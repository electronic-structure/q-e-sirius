set(src_elph_lib
	elph_bcast_input.f90  
	elphcom.f90  
	elph_compute_mel.f90  
	elph_load_q.f90  
	elph_prepare_q.f90  
	elph_q_points.f90  
	elph_readin.f90
	)

qe_add_library(qe_crpa_elph ${src_elph_lib})
target_link_libraries(qe_crpa_elph
    PRIVATE qe_crpa_w90
	    qe_phonon_ph
    	    qe_pw
            qe_lr_modules
            qe_modules
	    qe_lax
            qe_fftx
            qe_upflib
            qe_xclib)

###########################################################
# elph.x
###########################################################
set(src_elph_prog elph_main.f90)
qe_add_executable(qe_elph_exe ${src_elph_prog})
set_target_properties(qe_elph_exe PROPERTIES OUTPUT_NAME elph.x)
target_link_libraries(qe_elph_exe PRIVATE qe_crpa qe_crpa_w90 qe_crpa_elph qe_pw qe_lax qe_modules qe_xclib  qe_phonon_ph qe_lr_modules)
