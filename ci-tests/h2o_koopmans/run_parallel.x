#!/bin/bash

set -ex

pw.x -i h2o.scf.in -use_qe_scf
kcw.x -i h2o.kcw-wann2kcw.in
kcw.x -i h2o.kcw-screen.in

if [[ $SLURM_PROCID == 0 ]]; then
    python3 ../kcw_diff.py kcw.ref.yml kcw.yml
fi
