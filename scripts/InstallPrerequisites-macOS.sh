#!/bin/bash

# This script installs the prerequsites to build the Shasta
# static executable on macOS.
# These are necessary for building the executable,
# but not for running it. 
# The executable has no runtime dependencies other than
# standard system libraries that are always available.

# Cmake is not preinstalled on all systems.
brew install cmake

# Install the Boost libraries.
brew install boost

# The libpng and zlib libraries are needed for some
# functionality in the http server.
brew install libpng
# Assume zlib is always available.
# brew install zlib 

# Install Blas and Lapack.
# brew install openblas
# brew install lapack
# brew ls --verbose openblas
# brew ls --verbose lapack

# Install SeqAn 2.4.0
brew install brewsci/bio/seqan@2

# Build the spoa library (static library only).

# Create a temporary directory and cd to it.
tmpDirectoryName=$(mktemp -d /tmp/shasta-build.XXXX)
echo $tmpDirectoryName
cd $tmpDirectoryName

# Get the code.
curl -L https://github.com/rvaser/spoa/releases/download/3.4.0/spoa-v3.4.0.tar.gz \
    -o spoa-v3.4.0.tar.gz
tar -xvf spoa-v3.4.0.tar.gz

# Build the static library.
mkdir build-static
cd build-static
cmake ../spoa-v3.4.0 -DBUILD_SHARED_LIBS=OFF -Dspoa_optimize_for_native=OFF
make -j all
make install
cd 

# Remove the temporary directory.
rm -rf $tmpDirectoryName



# This builds and installs OGDF.
# Both the static and shared version of the library are built.
# It must be run under sudo.

# This installs the following:
# - Include file directories:
#       /usr/local/include/ogdf
#       /usr/local/include/coin
# - Static libraries:
#       /usr/local/lib/x86_64-linux-gnu/libOGDF.a
#       /usr/local/lib/x86_64-linux-gnu/libCOIN.a
# - Shared libraries:
#       /usr/local/lib/x86_64-linux-gnu/libOGDF.so
#       /usr/local/lib/x86_64-linux-gnu/libCOIN.so
# - Documentation:
#       /usr/local/share/doc/libogdf
# - Miscellaneous other files 
#       /usr/local/lib/x86_64-linux-gnu/cmake/OGDF

# Create a temporary directory to make the OGDF build.
tmpDirectoryName=$(mktemp --directory --tmpdir)
cd $tmpDirectoryName

# Get the code.
git clone https://github.com/ogdf/ogdf.git

# Build the shared library.
mkdir sharedBuild
cd sharedBuild
cmake ../ogdf -DBUILD_SHARED_LIBS=ON
make -j 8 all
make install
cd ..

# Build the static library.
mkdir staticBuild
cd staticBuild
cmake ../ogdf
make -j 8 all
make install

# Remove the temporary directory.
rm -rf $tmpDirectoryName



echo "==================================="
echo " Relevant software versions"
echo "==================================="
echo "$(brew list --versions cmake boost libpng seqan@2)"
echo "$(c++ --version | head -n 1)"
echo "==================================="


