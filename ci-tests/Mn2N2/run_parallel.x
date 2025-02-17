#!/bin/bash

set -ex

pw.x -i scf.in -use_qe_scf -npool 4
hp.x -i hp.in -npool 4

if [[ $SLURM_PROCID == 0 ]]; then
    cat Mn2N2.Hubbard_parameters.dat
    python3 ../hp_diff.py hp.ref.yml hp.yml
fi

