include(FindPackageHandleStandardArgs)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(_sciname_pe "gnu")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(_sciname_pe "intel")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set(_sciname_pe "cray")
else()
  message(${CMAKE_Fortran_COMPILER_ID})
  message(FATAL_ERROR "Unknown compiler. When using libsci use either GNU or INTEL compiler")
endif()

if (QE_ENABLE_MPI)
  if (QE_ENABLE_OPEMP)
    set(_sciname "sci_${_sciname_pe}_mpi_mp")
  else()
    set(_sciname "sci_${_sciname_pe}_mpi")
  endif()
else()
  # no MPI
  if (QE_ENABLE_OPEMP)
    set(_sciname "sci_${_sciname_pe}_mp")
  else()
    set(_sciname "sci_${_sciname_pe}")
  endif()
endif()

find_library(CRAY_LIBSCI_LIBRARIES SHARED
  NAMES ${_sciname}
  HINTS
  ${_SCALAPACK_LIBRARY_DIRS}
  ENV SCALAPACKROOT
  ENV CRAY_LIBSCI_PREFIX_DIR
  PATH_SUFFIXES lib
  DOC "scalapack library path")

find_package_handle_standard_args(CRAY_LIBSCI DEFAULT_MSG CRAY_LIBSCI_LIBRARIES)
