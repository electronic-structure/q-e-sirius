#!/bin/bash

set -ex

pw.x -i LiCoO2.scf.in -use_qe_scf -npool 2
hp.x -i LiCoO2.hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat LiCoO2.Hubbard_parameters.dat
    python3 ../hp_diff.py hp.ref.yml hp.yml
fi

