#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/Ni $PWD/Ni
else
    sleep 10
fi

cd $PWD/Ni
/apps/bin/pw.x -i Ni.scf.in -use_qe_scf -npool 3
/apps/bin/hp.x -i Ni.hp.in -npool 3
