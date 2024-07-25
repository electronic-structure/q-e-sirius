#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/CaFeO2 $PWD/CaFeO2
else
    sleep 10
fi

cd $PWD/CaFeO2
/apps/bin/pw.x -i CaFeO2.scf1.in -use_qe_scf -npool 6
/apps/bin/pw.x -i CaFeO2.scf2.in -use_qe_scf -npool 6
/apps/bin/hp.x -i CaFeO2.hp.in -npool 6

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/CaFeO2/CaFeO2.Hubbard_parameters.dat
    python3 /qe-src/ci-tests/hp_diff.py /qe-src/ci-tests/CaFeO2/hp.ref.yml $PWD/CaFeO2/hp.yml
else
    sleep 10
fi

