#!/bin/bash

set -xeuo pipefail

export SPACK_SYSTEM_CONFIG_PATH=/user-environment/config

SPACK_INSTALL_TREE=/dev/shm/spack-install

# make sure we keep the stage direcorty

spack env create -d ./spack-env
# add local repository with current sirius recipe
spack -e ./spack-env repo add $REPO
spack -e ./spack-env config add "packages:all:variants:[cuda_arch=${CUDA_ARCH},+cuda]"
# debug
cat ./spack-env/spack.yaml

# workaround, first command fails asking to update config format, doesn't make any sense, cannot reproduce on cli
#spack -e ./spack-env config add config:install_tree:$SPACK_INSTALL_TREE
yq w -i ./spack-env/spack.yaml 'spack.config.install_tree' $SPACK_INSTALL_TREE

spack -e ./spack-env add $SPEC


# build sirius from source
spack -e ./spack-env develop -p $PWD q-e-sirius@=develop-ristretto ^sirius@develop+cuda

# display spack.yaml
cat ./spack-env/spack.yaml

spack -e ./spack-env concretize
spack -e ./spack-env install

# create a symlink to spack build directory (keep in artifacts)
tar -cf installdir.tar $SPACK_INSTALL_TREE
