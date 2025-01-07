ARG BASE_IMAGE
FROM $BASE_IMAGE

# copy source files of the pull request into container
COPY . /qe-src

RUN spack -e /build-env-gcc add "q-e-sirius@develop-ristretto +libxc ^fftw+openmp ^sirius" && \
    spack -e /build-env-gcc develop -p /qe-src q-e-sirius@develop-ristretto && \
    spack -e /build-env-gcc concretize -f && \
    spack -e /build-env-gcc install

RUN spack clean --all

RUN spack arch

ENV PATH="$PATH:/apps/bin"

#RUN apt-get install -y pip && pip install reframe-hpc
#
#RUN git clone https://github.com/electronic-structure/qe-verification-tests.git
#ENV OMP_NUM_THREADS=1
#RUN cd qe-verification-tests/verification && reframe -c ./checks/ -r --system=localhost -C ./checks/config.py --tag serial

