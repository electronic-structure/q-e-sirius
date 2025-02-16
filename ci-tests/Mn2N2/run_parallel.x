#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r ./ci-tests/Mn2N2 $PWD/Mn2N2
else
    sleep 10
fi

cd $PWD/Mn2N2
pw.x -i scf.in -use_qe_scf -npool 5
hp.x -i hp.in -npool 5

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/Mn2N2.Hubbard_parameters.dat
    python3 ./ci-tests/hp_diff.py ./ci-tests/Mn2N2/hp.ref.yml $PWD/hp.yml
fi

