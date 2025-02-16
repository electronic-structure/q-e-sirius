#!/bin/bash

set -ex
cd $PWD/CaFeO2
pw.x -i CaFeO2.scf1.in -use_qe_scf -npool 2
pw.x -i CaFeO2.scf2.in -npool 2
hp.x -i CaFeO2.hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat CaFeO2.Hubbard_parameters.dat
    python3 ../hp_diff.py hp.ref.yml hp.yml
fi

