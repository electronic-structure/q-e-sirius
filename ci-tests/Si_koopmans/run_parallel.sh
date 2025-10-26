#!/bin/bash

set -ex

#scf
srun ${QE_PATH}/pw.x -npool 2 -i Si.scf.in -use_qe_scf

#wannier pp
srun -n1    ${QE_PATH}/wannier90.x -pp Si
srun -n1    ${QE_PATH}/wannier90.x -pp Si_emp

#pw2wannier
srun -n1 ${QE_PATH}/pw2wannier90.x -i Si.pw2wann.in
srun -n1 ${QE_PATH}/pw2wannier90.x -i Si_emp.pw2wann.in

#wannier
srun -n1 ${QE_PATH}/wannier90.x Si
srun -n1 cat Si.wout

srun -n1 ${QE_PATH}/wannier90.x Si_emp
srun -n1 cat Si_emp.wout

#kcw
srun ${QE_PATH}/kcw.x -i Si.kcw-wann2kcw.in
srun ${QE_PATH}/kcw.x -npool 2 -i Si.kcw-screen.in
source /user-environment/venv/bin/activate
srun -n1    python3 ../kcw_diff.py kcw.ref.yml kcw.yml
