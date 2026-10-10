set(src_w90_lib
	w90_interface.f90    
	w90_multiply_iqr.f90      
	w90_real_space.f90    
	w90_wan_Rr_buffer.f90  
	w90_interpolate.f90  
	w90_readwrite_wan_Rr.f90  
	w90_shift_center.f90  
	w90_wigner.f90
	)

qe_add_library(qe_crpa_w90 ${src_w90_lib})
target_link_libraries(qe_crpa_w90
    PRIVATE qe_pw
            qe_lr_modules
            qe_modules
	    qe_lax
            qe_fftx
            qe_upflib
            qe_xclib)

###########################################################
# wannier2pw.x
###########################################################
set(src_w90_prog wannier2pw.f90)
qe_add_executable(qe_wannier2pw_exe ${src_w90_prog})
set_target_properties(qe_wannier2pw_exe PROPERTIES OUTPUT_NAME wannier2pw.x)

target_include_directories(qe_wannier2pw_exe PRIVATE
    ${CMAKE_SOURCE_DIR}/LAXlib
)

target_link_libraries(qe_wannier2pw_exe PRIVATE qe_crpa_w90 qe_pw qe_lax qe_modules qe_xclib)
