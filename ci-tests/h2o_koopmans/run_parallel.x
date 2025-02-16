#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r ./ci-tests/h2o_koopmans $PWD/h2o_koopmans_parallel
else
    sleep 10
fi
cd $PWD/h2o_koopmans_parallel
pw.x -i h2o.scf.in -use_qe_scf
kcw.x -i h2o.kcw-wann2kcw.in
kcw.x -i h2o.kcw-screen.in 

if [[ $SLURM_PROCID == 0 ]]; then
    python3 ./ci-tests/kcw_diff.py ./ci-tests/h2o_koopmans/kcw.ref.yml $PWD/kcw.yml
fi
