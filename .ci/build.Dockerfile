FROM docker.io/electronicstructure/sirius

RUN spack external find --all --scope system --not-buildable

## show the spack's spec
#RUN spack spec -I $SPECDEV
#
#RUN spack env create --with-view /opt/sirius sirius-env
#RUN spack -e sirius-env add $SPECDEV
#
## copy source files of the pull request into container
COPY . /qe-src

ENV SPEC_QE="q-e-sirius@develop-ristretto ^fftw+openmp ^${SPEC_GCC_CPU}"

RUN spack spec -I $SPEC_QE

RUN spack install --only=dependencies $SPEC_QE &> /dev/null

RUN spack --color always dev-build --source-path /qe-src $SPEC_QE
#
## we need a fixed name for the build directory
## here is a hacky workaround to link ./spack-build-{hash} to ./spack-build
#RUN cd /sirius-src && ln -s $(find . -name "spack-build-*" -type d) spack-build
