FROM docker.io/electronicstructure/sirius

RUN spack spec -I $SPEC_GCC_CPU

RUN spack install $SPEC_GCC_CPU


## show the spack's spec
#RUN spack spec -I $SPECDEV
#
#RUN spack env create --with-view /opt/sirius sirius-env
#RUN spack -e sirius-env add $SPECDEV
#
## copy source files of the pull request into container
COPY . /qe-src

RUN spack --color always dev-build --source-path /qe-src q-e-sirius@develop-ristretto
#
## we need a fixed name for the build directory
## here is a hacky workaround to link ./spack-build-{hash} to ./spack-build
#RUN cd /sirius-src && ln -s $(find . -name "spack-build-*" -type d) spack-build
