###########################################################
# External FFT library that FFTXLIB depends on
# The following targets will be defined:
add_library(qe_ext_fft INTERFACE)
qe_install_targets(qe_ext_fft)
###########################################################
if(NOT QE_FFTW_INTERNAL)
    find_package(VendorPerfLibs COMPONENTS fft)
    if(NOT VendorPerfLibs_FOUND)
        message(FATAL_ERROR "Failed to find an external FFTW library. "
                            "Alternatively, '-DQE_FFTW_INTERNAL=ON' may be used to enable reference FFTW "
                            "at a performance loss compared to optimized libraries.")
    endif()

    target_link_libraries(qe_ext_fft INTERFACE VPL::fft)

        if(VPL_fft_ID STREQUAL "IntelMKL")
            qe_add_global_compile_definitions(__DFTI)
            set(qe_ext_fft_wrappers fft_scalar.DFTI.f90)
        else()
            qe_add_global_compile_definitions(__FFTW3)
            set(qe_ext_fft_wrappers fft_scalar.FFTW3.f90)
        endif()
else()
    message(STATUS "QE internal implementation of FFTW (FFTXLib)")
    qe_add_global_compile_definitions(__FFTW)
    set(qe_ext_fft_wrappers fft_scalar.FFTW.f90)
endif()

add_subdirectory(src)

###########################################################
# Tests
# TODO move all tests to a proper location
###########################################################
if(QE_ENABLE_TEST)
    add_subdirectory(tests)
endif()
