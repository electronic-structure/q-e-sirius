#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r ./ci-tests/FeO $PWD/FeO
else
    sleep 10
fi

cd $PWD/FeO
pw.x -i scf.in -use_qe_scf -npool 2
hp.x -i hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/FeO.Hubbard_parameters.dat
    python3 ../hp_diff.py hp.ref.yml hp.yml
fi

