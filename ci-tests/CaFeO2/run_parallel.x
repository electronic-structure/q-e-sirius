#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r ./ci-tests/CaFeO2 $PWD/CaFeO2
else
    sleep 10
fi

cd $PWD/CaFeO2
pw.x -i CaFeO2.scf1.in -use_qe_scf -npool 2
pw.x -i CaFeO2.scf2.in -npool 2
hp.x -i CaFeO2.hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/CaFeO2.Hubbard_parameters.dat
    python3 ./ci-tests/hp_diff.py ./ci-tests/CaFeO2/hp.ref.yml $PWD/hp.yml
fi

