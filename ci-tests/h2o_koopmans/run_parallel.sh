#!/bin/bash

set -ex

srun ${QE_PATH}/pw.x -i h2o.scf.in -use_qe_scf
srun ${QE_PATH}/kcw.x -i h2o.kcw-wann2kcw.in
srun ${QE_PATH}/kcw.x -i h2o.kcw-screen.in

source /user-environment/venv/bin/activate
srun -n1 python3 ../kcw_diff.py kcw.ref.yml kcw.yml
