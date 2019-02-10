#!/bin/bash

# Copyright 2018 M. Riechert and D. Meyer. Licensed under the MIT License.

set -ex

if [[ $OS_NAME == 'Linux' ]]; then

    sudo apt-get update
    sudo apt-get install gfortran libnetcdf-dev libnetcdff-dev netcdf-bin libpng-dev libjasper-dev

    if [ $BUILD_SYSTEM == 'Make' ]; then
        sudo apt-get install csh m4 libhdf5-serial-dev
    fi

    if [[ $MODE == dm* ]]; then
        sudo apt-get install libmpich-dev
    fi

    nc-config --all
    nf-config --all

elif [[ $OS_NAME == 'macOS' ]]; then

    # If c++ is already present on the system, homebrew fails to install gcc.
    # Use the `-f` flag in case c++ is not present to avoid errors.
    rm -f /usr/local/include/c++

    # disable automatic cleanup, just takes time
    export HOMEBREW_NO_INSTALL_CLEANUP=1

    brew update
    brew install coreutils gcc netcdf jasper libpng

    if [[ $MODE == dm* ]]; then
        brew install mpich
    fi

    nc-config --all

    # Homebrew installs the CMake version of netcdf which doesn't have nf-config support:
    # "nf-config not yet implemented for cmake builds".
    # This means WRF-Make won't enable NetCDF v4 support. For some reason, symlinking nc-config
    # to nf-config (as done for Ubuntu, see above) doesn't work here:
    # "/usr/local/bin/nf-config: fork: Resource temporarily unavailable"
    which nf-config
    #nf-config --has-nc4

else
    echo "The environment is not recognised"
fi

if [[ $MODE == dm* ]]; then
    mpif90 -v
fi
