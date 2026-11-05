#!/bin/bash

export OMP_NUM_THREADS=2

set -ex

#scf
pw.x -npool 1 -i Si.scf.in -use_qe_scf

#wannier pp
wannier90.x -pp Si

#pw2wannier
pw2wannier90.x -i Si.pw2wann.in

#wannier
wannier90.x Si

#crpa
crpa.x -i Si.crpa.in

python3 ../crpa_diff.py crpa_iq1.ref.yml crpa_iq1.yml
