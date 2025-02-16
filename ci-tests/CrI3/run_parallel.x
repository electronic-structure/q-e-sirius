#!/bin/bash

set -ex

pw.x -i CrI3.scf1.in -use_qe_scf -npool 2
pw.x -i CrI3.scf2.in -npool 2
hp.x -i CrI3.hp.in -npool 2
