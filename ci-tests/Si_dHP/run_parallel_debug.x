#!/bin/bash

export OMP_NUM_THREADS=2

set -ex

#scf
srun -N1 -n2 -c2 pw.x -npool 2 -i Si.scf.in -use_qe_scf

#wannier pp
srun -N1 -n1 -c2   wannier90.x -pp Si

#pw2wannier
srun -N1 -n1 -c2 pw2wannier90.x -i Si.pw2wann.in

#wannier
srun -N1 -n1 -c2 wannier90.x Si
srun -N1 -n1 -c2 cat Si.wout

#crpa
srun -N1 -n1 -c2 crpa.x -i Si.crpa.in
#source /user-environment/venv/bin/activate
#srun -n1    python3 ../crpa_diff.py crpa_iq1.ref.yml crpa_iq1.yml
