#!/bin/bash

# This installs all packages needed to build Shasta or to run
# the dynamic executable. The static executable has no prerequisites.

apt-get update
apt install -y git
apt install -y g++
apt install -y cmake
apt install -y libboost-all-dev
apt install -y libpng-dev
apt install -y graphviz
apt install -y gnuplot
apt install -y ncbi-blast+
apt install -y python3
apt install -y python3-pip
pip3 install pybind11



# SeqAn 2.4.0 is not available in the Ubuntu repository.
# Download it from GitHub, then install it.
apt install -y curl
curl -L https://github.com/seqan/seqan/releases/download/seqan-v2.4.0/seqan-library-2.4.0.deb \
    -o /dev/shm/seqan-library-2.4.0.deb
apt install /dev/shm/seqan-library-2.4.0.deb
rm /dev/shm/seqan-library-2.4.0.deb



# The spoa library is not available in the stable Ubuntu repository yet.
# Download it from GitHub, then install it.

# Create a temporary directory.
tmpDirectoryName=$(mktemp --directory --tmpdir)
echo $tmpDirectoryName
cd $tmpDirectoryName

# Get the code.
curl -L https://github.com/rvaser/spoa/releases/download/3.4.0/spoa-v3.4.0.tar.gz \
    -o spoa-v3.4.0.tar.gz
tar -xvf spoa-v3.4.0.tar.gz

ubuntuVersion=$(cat /etc/os-release | grep VERSION_ID | cut -f2 -d'=')
spoaGenDispatchFlag="ON"
if [ $ubuntuVersion \< "\"18.04\"" ] ; then
    # SPOA's cpu dispatching code is tested on newer Linux versions.
    spoaGenDispatchFlag="OFF"
fi

# Build the shared library.
mkdir build
cd build
cmake ../spoa-v3.4.0 -DBUILD_SHARED_LIBS=ON -Dspoa_generate_dispatch=$spoaGenDispatchFlag
make -j all
make install

# Build the static library.
cd ..
mkdir build-static
cd build-static
cmake ../spoa-v3.4.0 -DBUILD_SHARED_LIBS=OFF -Dspoa_generate_dispatch=$spoaGenDispatchFlag
make -j all
make install
cd 
rm -rf $tmpDirectoryName


# Make sure the newly created libraries are immediately visible to the loader.
ldconfig


