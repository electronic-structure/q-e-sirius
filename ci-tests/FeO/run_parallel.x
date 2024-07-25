#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/FeO $PWD/FeO
else
    sleep 10
fi

cd $PWD/FeO
/apps/bin/pw.x -i scf.in -use_qe_scf -npool 2
/apps/bin/hp.x -i hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/FeO.Hubbard_parameters.dat
    python3 /qe-src/ci-tests/hp_diff.py /qe-src/ci-tests/FeO/hp.ref.yml $PWD/hp.yml
fi

