#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/Si_koopmans $PWD/Si_koopmans_parallel
else
    sleep 10
fi

cd $PWD/Si_koopmans_parallel
cd $PWD/Si_koopmans
/apps/bin/pw.x -npool 2 -i Si.scf.in -use_qe_scf
/apps/bin/wannier90.x -pp Si
/apps/bin/pw2wannier90.x -i Si.pw2wann.in
/apps/bin/wannier90.x Si
cat Si.wout
/apps/bin/wannier90.x -pp Si_emp
/apps/bin/pw2wannier90.x -i Si_emp.pw2wann.in
/apps/bin/wannier90.x Si_emp
cat Si_emp.wout
/apps/bin/kcw.x -i Si.kcw-wann2kcw.in
/apps/bin/kcw.x -npool 2 -i Si.kcw-screen.in 

