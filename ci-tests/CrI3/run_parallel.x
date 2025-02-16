#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r ./ci-tests/CrI3 $PWD/CrI3
else
    sleep 10
fi

cd $PWD/CrI3
pw.x -i CrI3.scf1.in -use_qe_scf -npool 2
pw.x -i CrI3.scf2.in -npool 2
hp.x -i CrI3.hp.in -npool 2
