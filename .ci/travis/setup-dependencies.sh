#!/bin/bash

# Copyright 2018 M. Riechert and D. Meyer. Licensed under the MIT License.

set -ex

if [[ $TRAVIS_OS_NAME == 'linux' ]]; then

    sudo apt-get update

    # Travis CI is currently on Ubuntu 14.04 which has netCDF 4.1.
    # All versions of netCDF <= 4.1 contain all components (incl. Fortran libraries),
    # whereas netCDF > 4.1 is split up into separate libraries.
    # As soon as Travis CI switches to a new Ubuntu version, the below has to be adapted.
    sudo apt-get install gfortran libnetcdf-dev libpng-dev libjasper-dev

    if [ $SYSTEM == "make" ]; then
        sudo apt-get install csh m4 libhdf5-serial-dev
    fi

    if [[ $MODE == dm* ]]; then
        sudo apt-get install libmpich-dev
    fi

elif [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # Homebrew fails during "brew link" due to c++ being already present.
    rm '/usr/local/include/c++'

    brew update
    brew install gcc netcdf jasper

    if [[ $MODE == dm* ]]; then
        brew install mpich
    fi

else
    echo "The environment is not recognised"
fi