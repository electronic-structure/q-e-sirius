#!/bin/bash

set -ex

srun -N1 -n4 -c2 pw.x -i h2o.scf.in -use_qe_scf
srun -N1 -n4 -c2 kcw.x -i h2o.kcw-wann2kcw.in
srun -N1 -n4 -c2 kcw.x -i h2o.kcw-screen.in

