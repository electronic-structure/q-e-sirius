###########################################################
# WANNIER90
###########################################################
add_library(qe_wannier90 INTERFACE)
qe_install_targets(qe_wannier90)
if(QE_WANNIER90_INTERNAL)
    message(STATUS "Installing Wannier90 via submodule")

    # Wannier90 v4's own CMake asks for a newer CMake than QE does, and finds
    # BLAS/LAPACK for itself with find_package(... REQUIRED). Both only bite on
    # this path, so check them here rather than letting configuration fail
    # somewhere inside the submodule.
    if(CMAKE_VERSION VERSION_LESS 3.25)
        message(FATAL_ERROR
            "Wannier90 v4 requires CMake >= 3.25 (QE itself requires >= 3.20). "
            "Upgrade CMake, or build against an external Wannier90 with "
            "-DQE_WANNIER90_INTERNAL=OFF -DWANNIER90_ROOT=<path>.")
    endif()

    qe_git_submodule_update(external/wannier90)

    if(BLA_VENDOR STREQUAL "NVPL")
        # Wannier90 v4's own CMakeLists.txt finds BLAS/LAPACK for itself with
        # find_package(BLAS/LAPACK REQUIRED), which know nothing about NVPL as
        # a BLA_VENDOR and would otherwise fail or silently link a different
        # BLAS/LAPACK than the rest of QE. Shadow those two stock Find modules,
        # for this add_subdirectory() only, with ones that hand Wannier90 the
        # same NVPL targets QE itself resolved above. Saved/restored around
        # add_subdirectory(wannier90) so the shim doesn't leak into the other
        # external/ plugins (mbd, d3q, pw2qmcpack, qe-gipaw) processed after it.
        set(_qe_wannier90_saved_module_path ${CMAKE_MODULE_PATH})
        list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake/nvpl_wannier90")
    endif()

    set(WANNIER90_SHARED_LIBS ${BUILD_SHARED_LIBS})
    # Wannier90's own install rules must run so that Wannier90_lib/Wannier90_post
    # end up in an export set; otherwise qe_wannier90 (which INTERFACE-links them)
    # cannot be exported as part of qeTargets.
    set(WANNIER90_INSTALL ON)
    set(WANNIER90_TEST OFF)
    # Wannier90 defaults WANNIER90_MPI to OFF, so it has to be told to follow QE:
    # a serial library in an MPI build would take the per-rank m_local array that
    # EPW hands it via w90_set_m_local for the whole matrix, and would not report
    # anything, since valid_communicator() is unconditionally true without MPI.
    # WANNIER90_MPIH selects mpif.h over "use mpi", mirroring the COMMS=mpih vs
    # mpi90 choice in install/make_wannier90.inc.in. Setting these as normal
    # variables wins over Wannier90's option() calls because its own
    # cmake_minimum_required leaves CMP0077 at NEW.
    set(WANNIER90_MPI ${QE_ENABLE_MPI})
    if(QE_ENABLE_MPI AND NOT QE_ENABLE_MPI_MODULE)
        set(WANNIER90_MPIH ON)
    endif()
    add_subdirectory(wannier90)

    # nvfortran 25.x aborts with "gen_llvm_expr(): no incoming ili" on the
    # recursive internal procedure in disentangle.F90 unless inlining is off.
    if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
        set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/wannier90/src/disentangle.F90
            TARGET_DIRECTORY Wannier90_lib
            PROPERTIES COMPILE_OPTIONS "-Mnoinline")
    endif()

    if(BLA_VENDOR STREQUAL "NVPL")
        set(CMAKE_MODULE_PATH ${_qe_wannier90_saved_module_path})
        unset(_qe_wannier90_saved_module_path)
    endif()

    target_link_libraries(qe_wannier90 INTERFACE Wannier90::wannier90)

    ###########################################################
    # w90chk2chk.x
    ###########################################################
    add_executable(qe_w90chk2chk_exe wannier90/src/w90chk2chk.F90)
    set_target_properties(qe_w90chk2chk_exe PROPERTIES OUTPUT_NAME w90chk2chk.x)
    target_link_libraries(qe_w90chk2chk_exe PRIVATE Wannier90::wannier90)

    ###########################################################

    add_custom_target(w90
        DEPENDS
            qe_wannier90 Wannier90_exe Wannier90_post qe_w90chk2chk_exe
        COMMENT
            "Maximally localised Wannier Functions")

    qe_install_targets(
        # Executables
        Wannier90_exe Wannier90_post qe_w90chk2chk_exe)
else()
    # Wannier90 v4 installs its own CMake package -- Wannier90Config.cmake with
    # the Wannier90::wannier90 target -- so the hand-written FindWannier90
    # module is no longer needed and the same target name serves both branches.
    # WANNIER90_ROOT is passed as a hint explicitly: CMake only honours the
    # upper-case spelling of <PackageName>_ROOT from 3.27 (CMP0144), and QE
    # still supports older CMake.
    find_package(Wannier90 CONFIG REQUIRED HINTS ${WANNIER90_ROOT})
    target_link_libraries(qe_wannier90 INTERFACE Wannier90::wannier90)
endif()
