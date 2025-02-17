#!/bin/bash

set -ex

#scf
srun pw.x -npool 2 -i scf.in -use_qe_scf

#nscf
srun pw.x -npool 2 -i nscf.in

#wannier pp
srun -n1 wannier90.x -pp wann
srun -n1 wannier90.x -pp wann_emp
#pw2wannier
srun pw2wannier90.x -i occ.pw2wann.in
srun pw2wannier90.x -i emp.pw2wann.in

#wannier
srun -n1 wannier90.x wann
srun -n1 cat wann.wout

srun -n1 wannier90.x wann_emp
srun -n1 cat wann_emp.wout

#kcw
srun kcw.x -i kcw-wann2kcw.in
srun kcw.x -npool 2 -i kcw-screen.in

# load python pyyaml
source /user-environment/venv/bin/activate

srun -n1 python3 ../kcw_diff.py kcw.ref.yml kcw.yml
