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
    python3 ./ci-tests/hp_diff.py ./ci-tests/FeO/hp.ref.yml $PWD/hp.yml
fi

