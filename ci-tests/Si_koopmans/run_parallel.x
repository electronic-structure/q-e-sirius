#!/bin/bash

set -ex

#scf
pw.x -npool 2 -i Si.scf.in -use_qe_scf

#wannier pp
if [[ $SLURM_PROCID == 0 ]]; then
    wannier90.x -pp Si
    wannier90.x -pp Si_emp
else
    sleep 20
fi

#pw2wannier
pw2wannier90.x -i Si.pw2wann.in
pw2wannier90.x -i Si_emp.pw2wann.in

#wannier
if [[ $SLURM_PROCID == 0 ]]; then
    wannier90.x Si
    cat Si.wout
else
    sleep 100
fi

if [[ $SLURM_PROCID == 0 ]]; then
    wannier90.x Si_emp
    cat Si_emp.wout
else
    sleep 100
fi

#kcw
kcw.x -i Si.kcw-wann2kcw.in
kcw.x -npool 2 -i Si.kcw-screen.in 

if [[ $SLURM_PROCID == 0 ]]; then
    python3 ../kcw_diff.py kcw.ref.yml kcw.yml
fi
