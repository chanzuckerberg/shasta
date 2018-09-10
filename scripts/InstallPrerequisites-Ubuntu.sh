#!/bin/sh

apt-get update
apt install -y git
apt install -y g++
apt install -y cmake
apt install -y libboost-all-dev
apt install -y libpng-dev
apt install -y graphviz
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
