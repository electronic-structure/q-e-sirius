#!/bin/bash

set -ex

pw.x -i NiO.scf1.in -npool 2
# pw.x -i NiO.scf2.in -use_qe_scf -npool 2
hp.x -i NiO.hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat NiO.Hubbard_parameters.dat
    python3 ../hp_diff.py hp.ref.yml hp.yml
fi

