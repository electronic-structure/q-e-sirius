#!/bin/bash

set -ex

srun pw.x -i NiO.scf1.in -npool 2
# pw.x -i NiO.scf2.in -use_qe_scf -npool 2
srun hp.x -i NiO.hp.in -npool 2

srun -n1 cat NiO.Hubbard_parameters.dat
source /user-environment/venv/bin/activate
srun -n1 python3 ../hp_diff.py hp.ref.yml hp.yml
