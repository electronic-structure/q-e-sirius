#!/bin/bash

set -ex

srun ${QE_PATH}/pw.x -i Ni.scf.in -npool 3
srun ${QE_PATH}/hp.x -i Ni.hp.in -npool 3

source /user-environment/venv/bin/activate
srun -n1 cat Ni.Hubbard_parameters.dat
srun -n1 python3 ../hp_diff.py hp.ref.yml hp.yml
