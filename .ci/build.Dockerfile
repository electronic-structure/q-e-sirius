ARG BASE_IMAGE
FROM $BASE_IMAGE

# copy source files of the pull request into container
COPY . /qe-src
#RUN git clone https://github.com/electronic-structure/q-e-sirius.git /qe-src

RUN spack -e /build-env-gcc add "q-e-sirius@develop-ristretto ^fftw+openmp ^sirius" && \
    spack -e /build-env-gcc develop -p /qe-src q-e-sirius@develop-ristretto && \
    spack -e /build-env-gcc concretize -f && \
    spack -e /build-env-gcc install

RUN spack clean --all

ENV PATH="$PATH:/apps/bin"

RUN apt-get install -y pip && pip install reframe-hpc

RUN git clone https://github.com/electronic-structure/qe-verification-tests.git
RUN cd qe-verification-tests/verification && reframe -c ./checks/ -r --system=localhost -C ./checks/config.py --tag serial

