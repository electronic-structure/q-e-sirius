#!/bin/bash

set -ex

echo $PATH

srun pw.x -i scf.in -use_qe_scf -npool 4
srun hp.x -i hp.in -npool 4

srun -n1 cat Mn2N2.Hubbard_parameters.dat
source /user-environment/venv/bin/activate
srun -n1 python3 ../hp_diff.py hp.ref.yml hp.yml
