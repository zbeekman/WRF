#!/bin/bash

# Copyright 2018 M. Riechert and D. Meyer. Licensed under the MIT License.

set -ex

SCRIPTDIR=$(dirname "$0")
cd $SCRIPTDIR/../..

mkdir build && cd build

# TODO remove -DMPI_* variables once these get auto-detected
# Note that -DMODE alone decides whether MPI is used or not.
CC=gcc FC=gfortran cmake -G "MSYS Makefiles" \
    -DCMAKE_INSTALL_PREFIX=install -DCMAKE_BUILD_TYPE=${BUILD_TYPE} .. \
    -DMODE=$MODE -DNESTING=${NESTING} \
    -DENABLE_GRIB1=${GRIB1} -DENABLE_GRIB2=${GRIB2} \
    -DMPI_INCLUDE_PATH=$MINGW_PREFIX/include -DMPI_C_LIBRARY="$MSMPI_LIB64/msmpi.lib" \
    -DMPI_Fortran_LIBRARY="$MSMPI_LIB64/msmpifec.lib" ..

# It sometimes happens that the compiler runs out of memory due to parallel compilation.
# The construction below means "try with 2 cores, and if it fails, try again with 1 core".
cmake --build . --target install -- -j2 \
    || cmake --build . --target install \
    || cmake --build . --target install