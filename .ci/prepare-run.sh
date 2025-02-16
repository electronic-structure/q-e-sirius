#!/bin/bash

# spack install tree from the artifact to /dev/shm
# and sync all ranks

module load python

if [ $SLURM_LOCALID -eq 0 ]; then
    tar xf ./installdir.tar -C /
    touch /dev/shm/unpack_done_$CI_JOB_ID
    module load uv
    uv venv env
    source ./env/bin/activate
    uv pip install pyyaml
fi
# wait for tar unpack
while [ ! -f /dev/shm/unpack_done_$CI_JOB_ID ]; do
    sleep 0.2
done

source ./env/bin/activate

export PATH=/dev/shm/spack-install/q-e-sirius-develop-ristretto/bin/:$PATH
