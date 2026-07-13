set(sources
    w90_utils.f90
    w90_interface.f90
    w90_shift_center.f90
    w90_multiply_iqr.f90
    w90_wigner.f90
    w90_real_space.f90
    w90_readwrite_wan_Rr.f90
    w90_wan_Rr_buffer.f90
    w90_interpolate.f90
    w90_hamiltonian.f90)
qe_enable_cuda_fortran("${sources}")
qe_add_library(qe_w90 ${sources})
target_link_libraries(qe_w90
    PRIVATE
        qe_openacc_fortran
        qe_openmp_fortran
        qe_mpi_fortran
        qe_modules
        qe_pw
        qe_upflib
        qe_fftx
        qe_lapack)

###########################################################
# wan2pw.x
# (renamed from the reference's wannier2pw.x to avoid colliding
# with the unrelated PP/src wannier2pw.x)
###########################################################
set(src_wan2pw_x wan2pw.f90)
qe_enable_cuda_fortran("${src_wan2pw_x}")
qe_add_executable(qe_w90_wan2pw_exe ${src_wan2pw_x})
set_target_properties(qe_w90_wan2pw_exe PROPERTIES OUTPUT_NAME wan2pw.x)
target_link_libraries(qe_w90_wan2pw_exe
    PRIVATE
        qe_w90
        qe_pw
        qe_modules
        qe_lax
        qe_fftx
        qe_upflib
        qe_xclib
        qe_lapack)

###########################################################

qe_install_targets(qe_w90 qe_w90_wan2pw_exe)
