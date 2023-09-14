#!/bin/bash
set -euo pipefail

#rpm -qi glibc-static libgfortran-static
cd $(realpath -- $(dirname -- "$BASH_SOURCE"))

ROOT=$PWD/software
mkdir -p $ROOT

# Build PETSc
export PETSC_DIR=$ROOT/petsc.git
export PETSC_ARCH=arch-linux-opt
if [ ! -d $PETSC_DIR ]; then
    git clone -b v3.19.5 --depth 1 https://gitlab.com/petsc/petsc.git $PETSC_DIR
fi
if [ ! -d $PETSC_DIR/$PETSC_ARCH ]; then
    OPTFLAGS="-g0 -O3 -mtune=generic"
    options=(
        --COPTFLAGS="${OPTFLAGS}"
        --FOPTFLAGS="${OPTFLAGS}"
        --CC_LINKER_FLAGS="-static"
        --FC_LINKER_FLAGS="-static"
        --LIBS="-lm"
        --with-cxx=0
        --with-mpi=0
        --with-x=0
        --with-debugging=0
        --with-shared-libraries=0
        --with-fortran-bindings=0
        --download-f2cblaslapack=1
    )
    pushd $PETSC_DIR
    ./configure ${options[@]}
    popd
fi
if [ ! -f $PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a ]; then
    make -C $PETSC_DIR
fi

# Build SSDC
export SSDC_DIR=$ROOT/ssdc.git
export SSDC_ARCH=$PETSC_ARCH
if [ ! -d $SSDC_DIR ]; then
    git clone --depth 1 git@github.com:aanslab/ssdc.git $SSDC_DIR
fi
if [ ! -d $SSDC_DIR/$PETSC_ARCH ]; then
    pushd $SSDC_DIR
    make config
    popd
fi
if [ ! -f $SSDC_DIR/$PETSC_ARCH/lib/libssdc.a ]; then
    make -C $SSDC_DIR
fi

# Build drivers
make
