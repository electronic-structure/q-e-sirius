#!/bin/bash

set -ex

pw.x -i Ni.scf.in -npool 3
hp.x -i Ni.hp.in -npool 3

if [[ $SLURM_PROCID == 0 ]]; then
    cat Ni.Hubbard_parameters.dat
    python3 ../hp_diff.py hp.ref.yml hp.yml
fi

