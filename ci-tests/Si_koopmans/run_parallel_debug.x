#!/bin/bash

export OMP_NUM_THREADS=2

set -ex

#scf
srun -N1 -n2 -c2 pw.x -npool 2 -i Si.scf.in -use_qe_scf

#wannier pp
srun -N1 -n1 -c2   wannier90.x -pp Si
srun -N1 -n1 -c2   wannier90.x -pp Si_emp

#pw2wannier
srun -N1 -n1 -c2 pw2wannier90.x -i Si.pw2wann.in
srun -N1 -n1 -c2 pw2wannier90.x -i Si_emp.pw2wann.in

#wannier
srun -N1 -n1 -c2 wannier90.x Si
srun -N1 -n1 -c2 cat Si.wout

srun -N1 -n1 -c2 wannier90.x Si_emp
srun -N1 -n1 -c2 cat Si_emp.wout

#kcw
srun -N1 -n1 -c2 kcw.x -i Si.kcw-wann2kcw.in
srun -N1 -n2 -c2 kcw.x -npool 2 -i Si.kcw-screen.in
#source /user-environment/venv/bin/activate
#srun -n1    python3 ../kcw_diff.py kcw.ref.yml kcw.yml
