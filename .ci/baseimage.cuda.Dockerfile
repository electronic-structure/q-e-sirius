FROM ubuntu:24.04 as builder

ARG CUDA_ARCH=90

ENV DEBIAN_FRONTEND=noninteractive \
    PATH="$PATH:/spack/bin"

ENV CMAKE_VERSION=3.31.2
ENV MPICH_VERSION=3.4.3

#
#ENV FORCE_UNSAFE_CONFIGURE 1
#

RUN apt-get -y update && apt-get install -y apt-utils

# install basic tools
RUN apt-get install -y gcc g++ gfortran git make unzip \
  vim wget pkg-config curl tcl m4 cpio automake autoconf apt-transport-https \
  ca-certificates gnupg software-properties-common patchelf meson

RUN apt-get -y upgrade

# install CMake
RUN wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-aarch64.tar.gz -O cmake.tar.gz && \
    tar zxvf cmake.tar.gz --strip-components=1 -C /usr

# get latest version of spack
RUN git clone https://github.com/spack/spack.git

# set the location of packages built by spack
RUN spack config --scope system add config:install_tree:root:/opt/local
# set cuda_arch for all packages
RUN spack config --scope system add packages:all:variants:cuda_arch=${CUDA_ARCH}

# find gcc and clang compilers
RUN spack compiler find --scope system
RUN spack external find --all --scope system --not-buildable bash perl sed gcc  m4 \
    tar xz bzip2 cpio \
    cmake gmake make ninja meson autoconf automake \
    binutils findutils diffutils coreutils git curl openssh openssl ncurses

RUN spack install mpich@${MPICH_VERSION}

# for the MPI hook
RUN echo $(spack find --format='{prefix.lib}' mpich) > /etc/ld.so.conf.d/mpich.conf
RUN ldconfig

# install openblas
RUN spack install openblas threads=openmp %gcc +fortran

# install libvdwxc
RUN spack install libvdwxc %gcc +mpi ^mpich@${MPICH_VERSION}

RUN spack install nlcglib@develop %gcc +cuda

# create environments for several configurations and install dependencies
RUN spack env create -d /build-env-gcc --with-view /apps && \
    spack -e /build-env-gcc add "sirius@develop %gcc build_type=RelWithDebInfo +cuda +scalapack +vdwxc +fortran +tests +pugixml ^openblas%gcc ^libxc%gcc ^mpich%gcc ^umpire~device_alloc" && \
    spack -e /build-env-gcc concretize && \
    spack -e /build-env-gcc install --only=dependencies --fail-fast



