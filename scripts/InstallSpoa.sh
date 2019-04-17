#!/bin/bash


# The spoa library is not available in the Ubuntu repository.
# Download it from GitHub, then install it.

# Create a temporary directory.
tmpDirectoryName=$(mktemp --directory --tmpdir)
echo $tmpDirectoryName
cd $tmpDirectoryName

# Get the code.
curl -L https://github.com/rvaser/spoa/releases/download/3.0.0/spoa-v3.0.0.tar.gz \
    -o spoa-v3.0.0.tar.gz
tar -xvf spoa-v3.0.0.tar.gz

# Build the shared library.
mkdir build
cd build
cmake ../spoa-v3.0.0 -DBUILD_SHARED_LIBS=ON
make -j all
make install

# Build the static library.
cd ..
mkdir build-static
cd build-static
cmake ../spoa-v3.0.0 -DBUILD_SHARED_LIBS=OFF -Dspoa_optimize_for_native=OFF
make -j all
make install
cd 
rm -rf $tmpDirectoryName


