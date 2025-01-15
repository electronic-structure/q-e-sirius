#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/Mn2N2 $PWD/Mn2N2
else
    sleep 10
fi

cd $PWD/Mn2N2
/apps/bin/pw.x -i scf.in -use_qe_scf -npool 5
/apps/bin/hp.x -i hp.in -npool 5

if [[ $SLURM_PROCID == 0 ]]; then
    cat $PWD/Mn2N2.Hubbard_parameters.dat
    python3 /qe-src/ci-tests/hp_diff.py /qe-src/ci-tests/Mn2N2/hp.ref.yml $PWD/hp.yml
fi

