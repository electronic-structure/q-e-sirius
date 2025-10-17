#!/bin/bash

set -ex

srun ${QE_PATH}/pw.x -i scf.in -use_qe_scf -npool 2
srun ${QE_PATH}/hp.x -i hp.in -npool 2

srun -n1 cat NW.Hubbard_parameters.dat
source /user-environment/venv/bin/activate
srun -n1 python3 ../hp_diff.py hp.ref.yml hp.yml
