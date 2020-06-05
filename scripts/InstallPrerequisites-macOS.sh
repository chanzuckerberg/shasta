#!/bin/bash

# This script installs the prerequsites to build the Shasta
# static executable on macOS.
# These are necessary for building the executable,
# but not for running it. 
# The executable has no runtime dependencies other than
# standard system libraries that are always available.



# Install the Boost libraries.
brew install boost

# The libpng and zlib libraries are needed for some
# functionality in the http server.
brew install libpng
# Assume zlib is always available.
# brew install zlib 

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
cmake ../spoa-v3.4.0 -DBUILD_SHARED_LIBS=OFF
make -j all
make install
cd 

# Remove the temporary directory.
rm -rf $tmpDirectoryName

