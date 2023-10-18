FROM docker.io/electronicstructure/sirius

# copy source files of the pull request into container
COPY . /qe-src

RUN spack env create -d /build-env  --with-view /apps
RUN spack -e /build-env add q-e-sirius@develop-ristretto ^fftw+openmp ^${SPEC_GCC_GPU}
RUN spack -e /build-env develop -p /qe-src q-e-sirius@develop-ristretto

RUN spack -e /build-env concretize
RUN spack -e /build-env install

RUN spack clean --all

