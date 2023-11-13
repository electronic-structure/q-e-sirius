#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/CrI3 $PWD/CrI3
else
    sleep 10
fi

cd $PWD/CrI3
/apps/bin/pw.x -i CrI3.scf1.in -use_qe_scf -npool 2
/apps/bin/pw.x -i CrI3.scf2.in -use_qe_scf -npool 2
/apps/bin/hp.x -i CrI3.hp.in -npool 2
