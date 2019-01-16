#!/bin/bash

# This gets the marginPhase code from GitHub (polisher_shasta branch),
# builds the MarginCore library, and installs it in standard system locations.
# This code is also in InstallPrerequisites-Ubuntu.sh.

# Do everything in a temporary directory.
tmpDirectoryName=$(mktemp --directory --tmpdir)
pushd $tmpDirectoryName

# Get the code.
git clone https://github.com/benedictpaten/marginPhase.git
pushd marginPhase
git checkout 9da58634125452da362c840eee06e57d7fd4d48a
git submodule update --init
popd

# Build it.
mkdir build
pushd build
cmake ../marginPhase
make MarginCore

# Install it in standard system locations.
# This installs the following:
# - Shared library /usr/local/lib/libMarginCore.so,
#   required to build and run Shasta.
# - Header files in /usr/local/include/marginPhase
sudo make install

# Remove our temporary directory.
popd
popd 
rm -rf $tmpDirectoryName

# Make sure the newly created library is immediately visible to the loader.
sudo ldconfig


