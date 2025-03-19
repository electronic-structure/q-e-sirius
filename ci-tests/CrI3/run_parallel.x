#!/bin/bash

set -ex

srun pw.x -i CrI3.scf1.in -use_qe_scf -npool 2
srun pw.x -i CrI3.scf2.in -npool 2
srun hp.x -i CrI3.hp.in -npool 2

srun -n1 cat CrI3.Hubbard_parameters.dat
source /user-environment/venv/bin/activate
srun -n1 python3 ../hp_diff.py hp.ref.yml hp.yml
