#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/LiCoO2 $PWD/LiCoO2
else
    sleep 10
fi

cd $PWD/LiCoO2
/apps/bin/pw.x -i LiCoO2.scf.in -use_qe_scf -npool 2
/apps/bin/hp.x -i LiCoO2.hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/LiCoO2/LiCoO2.Hubbard_parameters.dat
    python3 /qe-src/ci-tests/hp_diff.py /qe-src/ci-tests/LiCoO2/hp.ref.yml $PWD/LiCoO2/hp.yml
else
    sleep 10
fi

