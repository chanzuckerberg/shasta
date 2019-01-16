#!/bin/bash



apt-get update
apt install -y git
apt install -y g++
apt install -y cmake
apt install -y libboost-all-dev
apt install -y libpng-dev
apt install -y graphviz
apt install -y ncbi-blast+
apt install -y hugepages
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



# The spoa library is not available in the Ubuntu repository.
# Download it from GitHub, then install it.
tmpDirectoryName=$(mktemp --directory --tmpdir)
echo $tmpDirectoryName
cd $tmpDirectoryName
curl -L https://github.com/rvaser/spoa/releases/download/1.1.3/spoa-v1.1.3.tar.gz \
    -o spoa-v1.1.3.tar.gz
tar -xvf spoa-v1.1.3.tar.gz
mkdir build
cd build
cmake ../spoa-v1.1.3 -DBUILD_SHARED_LIBS=ON
make -j all
make install
cd 
rm -rf $tmpDirectoryName



# This gets the marginPhase code from GitHub (polisher_shasta branch),
# builds the MarginCore library, and installs it in standard system locations.
# This code is also in InstallMarginCore.sh.

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



# Make sure the newly created libraries are immediately visible to the loader.
ldconfig


