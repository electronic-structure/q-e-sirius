#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/NiO $PWD/NiO
else
    sleep 10
fi

cd $PWD/NiO
/apps/bin/pw.x -i NiO.scf1.in -use_qe_scf -npool 4
/apps/bin/pw.x -i NiO.scf2.in -use_qe_scf -npool 4
/apps/bin/hp.x -i NiO.hp.in -npool 4

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/NiO.Hubbard_parameters.dat
    python3 /qe-src/ci-tests/hp_diff.py /qe-src/ci-tests/NiO/hp.ref.yml $PWD/hp.yml
fi

