#!/bin/bash

export OMP_NUM_THREADS=2

set -ex

#scf
pw.x -npool 1 -i Si.scf.in -use_qe_scf

#wannier pp
wannier90.x -pp Si
wannier90.x -pp Si_emp

#pw2wannier
pw2wannier90.x -i Si.pw2wann.in
pw2wannier90.x -i Si_emp.pw2wann.in

#wannier
wannier90.x Si
wannier90.x Si_emp

#kcw
kcw.x -i Si.kcw-wann2kcw.in
kcw.x -i Si.kcw-screen.in
