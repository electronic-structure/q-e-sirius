#!/bin/bash

set -ex

srun pw.x -i CrI3.scf1.in -use_qe_scf -npool 2
srun pw.x -i CrI3.scf2.in -npool 2
srun hp.x -i CrI3.hp.in -npool 2
