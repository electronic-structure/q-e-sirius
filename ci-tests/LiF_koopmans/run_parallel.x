#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/LiF_koopmans $PWD/LiF_koopmans_parallel
else
    sleep 10
fi

cd $PWD/LiF_koopmans_parallel

#scf
/apps/bin/pw.x -npool 2 -i scf.in -use_qe_scf

#nscf
/apps/bin/pw.x -npool 2 -i nscf.in

#wannier pp
if [[ $SLURM_PROCID == 0 ]]; then
    /apps/bin/wannier90.x -pp wann 
    /apps/bin/wannier90.x -pp wann_emp
else
    sleep 20
fi

#pw2wannier
/apps/bin/pw2wannier90.x -i occ.pw2wann.in
/apps/bin/pw2wannier90.x -i emp.pw2wann.in

#wannier
if [[ $SLURM_PROCID == 0 ]]; then
    /apps/bin/wannier90.x wann
    cat wann.wout
else
    sleep 100
fi

if [[ $SLURM_PROCID == 0 ]]; then
    /apps/bin/wannier90.x wann_emp
    cat wann_emp.wout
else
    sleep 100
fi

#kcw
/apps/bin/kcw.x -i kcw-wann2kcw.in
/apps/bin/kcw.x -npool 2 -i kcw-screen.in 

if [[ $SLURM_PROCID == 0 ]]; then
    python3 /qe-src/ci-tests/kcw_diff.py /qe-src/ci-tests/LiF_koopmans/kcw.ref.yml $PWD/kcw.yml
fi
