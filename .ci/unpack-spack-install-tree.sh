#!/bin/bash

# spack install tree from the artifact to /dev/shm
# and sync all ranks

# unpack if running within sbatch and also on baremetal (SLURM_LOCALID undef)
if [ -z "$SLURM_LOCALID" ] || [ $SLURM_LOCALID -eq 0 ]; then
    tar xf ./installdir.tar -C /
    touch /dev/shm/unpack_done_$CI_JOB_ID
fi
# wait for tar unpack
while [ ! -f /dev/shm/unpack_done_$CI_JOB_ID ]; do
    sleep 0.2
done

