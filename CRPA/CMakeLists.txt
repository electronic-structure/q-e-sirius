add_subdirectory(wannier)
add_subdirectory(src)
add_subdirectory(elph)

###########################################################

qe_install_targets(qe_crpa qe_crpa_exe qe_crpa_w90  qe_wannier2pw_exe  qe_elph_exe)

add_custom_target(
    crpa
    DEPENDS pw w90 qe_crpa_exe  qe_wannier2pw_exe  qe_elph_exe
    COMMENT "CRPA implementation")

