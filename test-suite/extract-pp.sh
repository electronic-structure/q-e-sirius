# Copyright (C) 2020 Quantum ESPRESSO Foundation
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
 
fname=$1

# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`

# PPACF
efock=`grep 'Fock energy' $fname | awk '{print $4}'`
exlda=`grep 'LDA Exchange' $fname | awk '{print $3}'`
eclda=`grep 'LDA Correlation' $fname | awk '{print $3}'`
exc=`grep 'Exchange + Correlation' $fname | awk '{print $4}'`
etcl=`grep 'T_c^LDA' $fname | awk '{print $2}'`
etnl=`grep 'T_c^nl' $fname | awk '{print $2}'`
ekc=`grep 'Kinetic-correlation Energy' $fname | awk '{print $3}'`
enl=`grep 'Non-local energy' $fname | awk '{print $4}'`

# PW2WANNIER90 in library mode: Wannier90 results, from the appended .wout
wf=`sed -n '/Final State/,/Sum of centres/p' $fname | grep 'WF centre and spread' \
    | sed 's/[(),]/ /g' | awk '{print $6, $7, $8, $9}'`
omegai=`grep 'Omega I  ' $fname | tail -1 | awk '{print $NF}'`
omegad=`grep 'Omega D  ' $fname | tail -1 | awk '{print $NF}'`
omegaod=`grep 'Omega OD ' $fname | tail -1 | awk '{print $NF}'`
omegatot=`grep 'Omega Total' $fname | tail -1 | awk '{print $NF}'`
disomegai=`grep 'Final Omega_I' $fname | tail -1 | awk '{print $3}'`

# Standalone Wannier90 chain: number of nearest neighbours from the .nnkp
# file, and the eigenvalues appended from the .eig file
nntot=`sed -n '/begin nnkpts/{n;p;}' $fname | awk '{print $1}'`
eig=`grep '^EIG ' $fname | awk '{print $2}'`

# PW2WANNIER90 with wan_mode='decompose': power spectra, from the appended .power files
pspec=`grep '^POWER ' $fname | awk '{print $2}'`
gcoeff=`grep '^GCOEFF ' $fname | awk '{print $2}'`

if test "$e1" != ""; then
        echo e1
        echo $e1
fi
if test "$efock" != ""; then
        echo efock
        echo $efock
fi
if test "$exlda" != ""; then
        echo exlda
        echo $exlda
fi
if test "$eclda" != ""; then
        echo eclda
        echo $eclda
fi
if test "$exc" != ""; then
        echo exc
        echo $exc
fi
if test "$etcl" != ""; then
        echo etcl
        echo $etcl
fi
if test "$etnl" != ""; then
        echo etnl
        echo $etnl
fi
if test "$ekc" != ""; then
        echo ekc
        echo $ekc
fi
if test "$enl" != ""; then
        echo enl
        echo $enl
fi
if test "$wf" != ""; then
        echo wfcx wfcy wfcz wfspread
        echo "$wf"
fi
if test "$omegai" != ""; then
        echo omegai
        echo $omegai
fi
if test "$omegad" != ""; then
        echo omegad
        echo $omegad
fi
if test "$omegaod" != ""; then
        echo omegaod
        echo $omegaod
fi
if test "$omegatot" != ""; then
        echo omegatot
        echo $omegatot
fi
if test "$disomegai" != ""; then
        echo disomegai
        echo $disomegai
fi
if test "$nntot" != ""; then
        echo nntot
        echo $nntot
fi
if test "$eig" != ""; then
        echo eig
        echo "$eig"
fi
if test "$pspec" != ""; then
        echo pspec
        echo "$pspec"
fi
if test "$gcoeff" != ""; then
        echo gcoeff
        echo "$gcoeff"
fi
