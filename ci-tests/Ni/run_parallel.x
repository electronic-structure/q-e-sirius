#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    spack arch
    cp -r ./ci-tests/Ni $PWD/Ni
else
    sleep 10
fi

cd $PWD/Ni
pw.x -i Ni.scf.in -npool 3
hp.x -i Ni.hp.in -npool 3

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/Ni.Hubbard_parameters.dat
    python3 ./ci-tests/hp_diff.py ./ci-tests/Ni/hp.ref.yml $PWD/hp.yml
fi

