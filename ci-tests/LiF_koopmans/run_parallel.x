#!/bin/bash

set -ex

#scf
pw.x -npool 2 -i scf.in -use_qe_scf

#nscf
pw.x -npool 2 -i nscf.in

#wannier pp
if [[ $SLURM_PROCID == 0 ]]; then
    wannier90.x -pp wann 
    wannier90.x -pp wann_emp
else
    sleep 20
fi

#pw2wannier
pw2wannier90.x -i occ.pw2wann.in
pw2wannier90.x -i emp.pw2wann.in

#wannier
if [[ $SLURM_PROCID == 0 ]]; then
    wannier90.x wann
    cat wann.wout
else
    sleep 100
fi

if [[ $SLURM_PROCID == 0 ]]; then
    wannier90.x wann_emp
    cat wann_emp.wout
else
    sleep 100
fi

#kcw
kcw.x -i kcw-wann2kcw.in
kcw.x -npool 2 -i kcw-screen.in

if [[ $SLURM_PROCID == 0 ]]; then
    python3 ../kcw_diff.py kcw.ref.yml kcw.yml
fi
