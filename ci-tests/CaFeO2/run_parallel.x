#!/bin/bash

set -ex
srun pw.x -i CaFeO2.scf1.in -use_qe_scf -npool 2
srun pw.x -i CaFeO2.scf2.in -npool 2
srun hp.x -i CaFeO2.hp.in -npool 2

srun -n1 cat CaFeO2.Hubbard_parameters.dat

source /user-environment/venv/bin/activate
srun -n1 python3 ../hp_diff.py hp.ref.yml hp.yml
