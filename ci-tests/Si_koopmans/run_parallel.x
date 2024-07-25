#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/Si_koopmans $PWD/Si_koopmans_parallel
else
    sleep 10
fi

cd $PWD/Si_koopmans_parallel
/apps/bin/pw.x -npool 2 -i Si.scf.in -use_qe_scf
if [[ $SLURM_PROCID == 0 ]]; then
    /apps/bin/wannier90.x -pp Si
    /apps/bin/wannier90.x -pp Si_emp
else
    sleep 20
fi

/apps/bin/pw2wannier90.x -i Si.pw2wann.in
/apps/bin/wannier90.x Si
if [[ $SLURM_PROCID == 0 ]]; then
    cat Si.wout
else
    sleep 20
fi
/apps/bin/pw2wannier90.x -i Si_emp.pw2wann.in
/apps/bin/wannier90.x Si_emp
if [[ $SLURM_PROCID == 0 ]]; then
    cat Si_emp.wout
else
    sleep 20
fi

/apps/bin/kcw.x -i Si.kcw-wann2kcw.in
/apps/bin/kcw.x -npool 2 -i Si.kcw-screen.in 

if [[ $SLURM_PROCID == 0 ]]; then
    python3 /qe-src/ci-tests/hp_diff.py /qe-src/ci-tests/Si_koopmans/kcw.ref.yml $PWD/Si_koopmans_parallel/kcw.yml
else
    sleep 10
fi
