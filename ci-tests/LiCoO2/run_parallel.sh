#!/bin/bash

set -ex

srun pw.x -i LiCoO2.scf.in -use_qe_scf -npool 2
srun hp.x -i LiCoO2.hp.in -npool 2

srun -n1 cat LiCoO2.Hubbard_parameters.dat

source /user-environment/venv/bin/activate
srun -n1 python3 ../hp_diff.py hp.ref.yml hp.yml
