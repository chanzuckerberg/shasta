#!/bin/bash

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

# Make sure the newly created libraries are immediately visible to the loader.
ldconfig


