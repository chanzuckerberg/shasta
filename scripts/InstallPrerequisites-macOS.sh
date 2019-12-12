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

# Build the spoa library (static library only).

# Create a temporary directory and cd to it.
tmpDirectoryName=$(mktemp -d /tmp/shasta-build.XXXX)
echo $tmpDirectoryName
cd $tmpDirectoryName

# Get the code.
curl -L https://github.com/rvaser/spoa/releases/download/3.0.0/spoa-v3.0.0.tar.gz \
    -o spoa-v3.0.0.tar.gz
tar -xvf spoa-v3.0.0.tar.gz

# Build the static library.
mkdir build-static
cd build-static
cmake ../spoa-v3.0.0 -DBUILD_SHARED_LIBS=OFF -Dspoa_optimize_for_native=OFF
make -j all
make install
cd ..  

# This gets the edlib code from GitHub,
# builds the edlib library, and installs it in standard system locations.

# Get the code.
git clone https://github.com/Martinsos/edlib.git
cd edlib
git checkout ba4272ba68fcdbe31cbc10853de1841701e4e60a 
cd ../

# Build it.
mkdir build
cd build
cmake -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true -DCMAKE_BUILD_TYPE=Release ../edlib
make edlib 

# Install it in standard system locations.
sudo make install

# Remove the temporary directory.
cd
rm -rf $tmpDirectoryName

