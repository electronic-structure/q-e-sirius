#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r ./ci-tests/NiO $PWD/NiO
else
    sleep 10
fi

cd $PWD/NiO
/apps/bin/pw.x -i NiO.scf1.in -npool 2
#/apps/bin/pw.x -i NiO.scf2.in -use_qe_scf -npool 2
/apps/bin/hp.x -i NiO.hp.in -npool 2

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/NiO.Hubbard_parameters.dat
    python3 ./ci-tests/hp_diff.py ./ci-tests/NiO/hp.ref.yml $PWD/hp.yml
fi

