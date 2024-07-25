#!/bin/bash

set -ex

if [[ $SLURM_PROCID == 0 ]]; then
    cp -r /qe-src/ci-tests/h2o_koopmans $PWD/h2o_koopmans_parallel
else
    sleep 10
fi
cd $PWD/h2o_koopmans_parallel
/apps/bin/pw.x -npool 2 -i h2o.scf.in -use_qe_scf
/apps/bin/kcw.x -i h2o.kcw-wann2kcw.in
/apps/bin/kcw.x -npool 2 -i h2o.kcw-screen.in 

if [[ $SLURM_PROCID == 0 ]]; then
    python3 /qe-src/ci-tests/hp_diff.py /qe-src/ci-tests/h2o_koopmans/kcw.ref.yml $PWD/kcw.yml
else
    sleep 10
fi
