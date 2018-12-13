#!/bin/sh



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
cd ..
rm -rf $tmpDirectoryName

# Make sure the newly created library is immediately visible to the loader.
ldconfig


