#!/bin/bash

set -ex

#scf
srun ${QE_PATH}/pw.x -npool 2 -i Si.scf.in -use_qe_scf

#wannier pp
srun -n1 ${QE_PATH}/wannier90.x -pp Si

#pw2wannier
srun -n1 ${QE_PATH}/pw2wannier90.x -i Si.pw2wann.in

#wannier
srun -n1 ${QE_PATH}/wannier90.x Si
srun -n1 cat Si.wout

#crpa
srun ${QE_PATH}/crpa.x -i Si.crpa.in
source /user-environment/venv/bin/activate
srun -n1    python3 ../crpa_diff.py crpa_iq1.ref.yml crpa_iq1.yml
