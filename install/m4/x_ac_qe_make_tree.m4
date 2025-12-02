# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_MAKE_TREE], [


AC_CONFIG_FILES([install/extlibs_makefile:install/extlibs_makefile])
AC_CONFIG_FILES([install/install_utils:install/install_utils])
AC_CONFIG_FILES([install/plugins_makefile:install/plugins_makefile])
AC_CONFIG_FILES([install/plugins_list:install/plugins_list])

AC_CONFIG_FILES([Makefile:Makefile])

AC_CONFIG_FILES([COUPLE/Makefile:COUPLE/Makefile])
AC_CONFIG_FILES([COUPLE/src/Makefile:COUPLE/src/Makefile])

AC_CONFIG_FILES([FFTXlib/Makefile:FFTXlib/Makefile])
AC_CONFIG_FILES([FFTXlib/tests/Makefile:FFTXlib/tests/Makefile])
AC_CONFIG_FILES([FFTXlib/src/Makefile:FFTXlib/src/Makefile])
AC_CONFIG_FILES([FFTXlib/src/make.depend:FFTXlib/src/make.depend])

AC_CONFIG_FILES([PIOUD/Makefile:PIOUD/Makefile])
AC_CONFIG_FILES([PIOUD/src/Makefile:PIOUD/src/Makefile])
AC_CONFIG_FILES([PIOUD/src/make.depend:PIOUD/src/make.depend])

AC_CONFIG_FILES([NEB/Makefile:NEB/Makefile])
AC_CONFIG_FILES([NEB/src/Makefile:NEB/src/Makefile])
AC_CONFIG_FILES([NEB/src/make.depend:NEB/src/make.depend])

AC_CONFIG_FILES([XSpectra/Makefile:XSpectra/Makefile])
AC_CONFIG_FILES([XSpectra/src/Makefile:XSpectra/src/Makefile])
AC_CONFIG_FILES([XSpectra/src/make.depend:XSpectra/src/make.depend])

AC_CONFIG_FILES([XClib/Makefile:XClib/Makefile])
AC_CONFIG_FILES([XClib/make.depend:XClib/make.depend])

AC_CONFIG_FILES([Modules/Makefile:Modules/Makefile])
AC_CONFIG_FILES([Modules/make.depend:Modules/make.depend])

AC_CONFIG_FILES([atomic/Makefile:atomic/Makefile])
AC_CONFIG_FILES([atomic/src/Makefile:atomic/src/Makefile])
AC_CONFIG_FILES([atomic/src/make.depend:atomic/src/make.depend])

AC_CONFIG_FILES([PP/Makefile:PP/Makefile])
AC_CONFIG_FILES([PP/src/Makefile:PP/src/Makefile])
AC_CONFIG_FILES([PP/src/make.depend:PP/src/make.depend])
AC_CONFIG_FILES([PP/simple_transport/src/Makefile:PP/simple_transport/src/Makefile])

AC_CONFIG_FILES([QEHeat/Makefile:QEHeat/Makefile])
AC_CONFIG_FILES([QEHeat/src/Makefile:QEHeat/src/Makefile])
AC_CONFIG_FILES([QEHeat/src/make.depend:QEHeat/src/make.depend])

AC_CONFIG_FILES([GWW/Makefile:GWW/Makefile])
AC_CONFIG_FILES([GWW/pw4gww/Makefile:GWW/pw4gww/Makefile])
AC_CONFIG_FILES([GWW/pw4gww/make.depend:GWW/pw4gww/make.depend])
AC_CONFIG_FILES([GWW/gww/Makefile:GWW/gww/Makefile])
AC_CONFIG_FILES([GWW/gww/make.depend:GWW/gww/make.depend])
AC_CONFIG_FILES([GWW/simple_bse/Makefile:GWW/simple_bse/Makefile])
AC_CONFIG_FILES([GWW/simple_bse/make.depend:GWW/simple_bse/make.depend])
AC_CONFIG_FILES([GWW/head/Makefile:GWW/head/Makefile])
AC_CONFIG_FILES([GWW/head/make.depend:GWW/head/make.depend])
AC_CONFIG_FILES([GWW/simple_ip/Makefile:GWW/simple_ip/Makefile])
AC_CONFIG_FILES([GWW/simple_ip/make.depend:GWW/simple_ip/make.depend])
AC_CONFIG_FILES([GWW/bse/Makefile:GWW/bse/Makefile])
AC_CONFIG_FILES([GWW/bse/make.depend:GWW/bse/make.depend])
AC_CONFIG_FILES([GWW/simple/Makefile:GWW/simple/Makefile])
AC_CONFIG_FILES([GWW/simple/make.depend:GWW/simple/make.depend])
AC_CONFIG_FILES([GWW/util/Makefile:GWW/util/Makefile])
AC_CONFIG_FILES([GWW/minpack/Makefile:GWW/minpack/Makefile])

AC_CONFIG_FILES([HP/Makefile:HP/Makefile])
AC_CONFIG_FILES([HP/src/Makefile:HP/src/Makefile])
AC_CONFIG_FILES([HP/src/make.depend:HP/src/make.depend])

AC_CONFIG_FILES([upflib/Makefile:upflib/Makefile])
AC_CONFIG_FILES([upflib/make.depend:upflib/make.depend])

AC_CONFIG_FILES([PW/Makefile:PW/Makefile])
AC_CONFIG_FILES([PW/tools/Makefile:PW/tools/Makefile])
AC_CONFIG_FILES([PW/tools/make.depend:PW/tools/make.depend])
AC_CONFIG_FILES([PW/src/Makefile:PW/src/Makefile])
AC_CONFIG_FILES([PW/src/make.depend:PW/src/make.depend])

AC_CONFIG_FILES([PWCOND/Makefile:PWCOND/Makefile])
AC_CONFIG_FILES([PWCOND/src/Makefile:PWCOND/src/Makefile])
AC_CONFIG_FILES([PWCOND/src/make.depend:PWCOND/src/make.depend])

AC_CONFIG_FILES([KS_Solvers/Makefile:KS_Solvers/Makefile])
AC_CONFIG_FILES([KS_Solvers/CG/Makefile:KS_Solvers/CG/Makefile])
AC_CONFIG_FILES([KS_Solvers/CG/make.depend:KS_Solvers/CG/make.depend])
AC_CONFIG_FILES([KS_Solvers/RMM/Makefile:KS_Solvers/RMM/Makefile])
AC_CONFIG_FILES([KS_Solvers/RMM/make.depend:KS_Solvers/RMM/make.depend])
AC_CONFIG_FILES([KS_Solvers/ParO/Makefile:KS_Solvers/ParO/Makefile])
AC_CONFIG_FILES([KS_Solvers/ParO/make.depend:KS_Solvers/ParO/make.depend])
AC_CONFIG_FILES([KS_Solvers/Davidson/Makefile:KS_Solvers/Davidson/Makefile])
AC_CONFIG_FILES([KS_Solvers/Davidson/make.depend:KS_Solvers/Davidson/make.depend])
AC_CONFIG_FILES([KS_Solvers/DENSE/Makefile:KS_Solvers/DENSE/Makefile])
AC_CONFIG_FILES([KS_Solvers/DENSE/make.depend:KS_Solvers/DENSE/make.depend])
AC_CONFIG_FILES([KS_Solvers/Davidson_RCI/Makefile:KS_Solvers/Davidson_RCI/Makefile])
AC_CONFIG_FILES([KS_Solvers/Davidson_RCI/make.depend:KS_Solvers/Davidson_RCI/make.depend])

AC_CONFIG_FILES([CPV/Makefile:CPV/Makefile])
AC_CONFIG_FILES([CPV/src/Makefile:CPV/src/Makefile])
AC_CONFIG_FILES([CPV/src/make.depend:CPV/src/make.depend])


AC_CONFIG_FILES([dft-d3/Makefile:dft-d3/Makefile])
AC_CONFIG_FILES([dft-d3/make.depend:dft-d3/make.depend])

# TG: TDDFPT/ColorCalculator/Install?
AC_CONFIG_FILES([TDDFPT/Makefile:TDDFPT/Makefile])
AC_CONFIG_FILES([TDDFPT/src/Makefile:TDDFPT/src/Makefile])
AC_CONFIG_FILES([TDDFPT/src/make.depend:TDDFPT/src/make.depend])

AC_CONFIG_FILES([UtilXlib/Makefile:UtilXlib/Makefile])
AC_CONFIG_FILES([UtilXlib/make.depend:UtilXlib/make.depend])
# TG: autotest.inc?
AC_CONFIG_FILES([UtilXlib/tests/Makefile:UtilXlib/tests/Makefile])

AC_CONFIG_FILES([LR_Modules/Makefile:LR_Modules/Makefile])
AC_CONFIG_FILES([LR_Modules/make.depend:LR_Modules/make.depend])

AC_CONFIG_FILES([LAXlib/Makefile:LAXlib/Makefile])
AC_CONFIG_FILES([LAXlib/make.depend:LAXlib/make.depend])
AC_CONFIG_FILES([LAXlib/tests/Makefile:LAXlib/tests/Makefile])

AC_CONFIG_FILES([PHonon/Makefile:PHonon/Makefile])
AC_CONFIG_FILES([PHonon/FD/Makefile:PHonon/FD/Makefile])
AC_CONFIG_FILES([PHonon/FD/make.depend:PHonon/FD/make.depend])
AC_CONFIG_FILES([PHonon/Gamma/Makefile:PHonon/Gamma/Makefile])
AC_CONFIG_FILES([PHonon/Gamma/make.depend:PHonon/Gamma/make.depend])
AC_CONFIG_FILES([PHonon/PH/Makefile:PHonon/PH/Makefile])
AC_CONFIG_FILES([PHonon/PH/make.depend:PHonon/PH/make.depend])

AC_CONFIG_FILES([KCW/Makefile:KCW/Makefile])
AC_CONFIG_FILES([KCW/PP/Makefile:KCW/PP/Makefile])
AC_CONFIG_FILES([KCW/PP/make.depend:KCW/PP/make.depend])
AC_CONFIG_FILES([KCW/src/Makefile:KCW/src/Makefile])
AC_CONFIG_FILES([KCW/src/make.depend:KCW/src/make.depend])

AC_CONFIG_FILES([EPW/Makefile:EPW/Makefile])
AC_CONFIG_FILES([EPW/irobjs/Makefile:EPW/irobjs/Makefile])
AC_CONFIG_FILES([EPW/ZG/src/Makefile:EPW/ZG/src/Makefile])
AC_CONFIG_FILES([EPW/ZG/src/make.depend:EPW/ZG/src/make.depend])
AC_CONFIG_FILES([EPW/src/Makefile:EPW/src/Makefile])
AC_CONFIG_FILES([EPW/src/make.depend:EPW/src/make.depend])

AC_CONFIG_FILES([archive/README.md:archive/README.md])

]
)
