#!/bin/bash

set -ex

# Copyright 2018 M. Riechert and D. Meyer. Licensed under the MIT License.

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

    # This is really Python 2.
    brew uninstall python || true

    # Install latest Python 3.
    brew update
    brew install python

    python3 -V
    pip3 -V

fi
