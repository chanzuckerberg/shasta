#!/bin/bash

# This installs all packages needed to build Shasta or to run
# the dynamic executable. The static executable has no prerequisites.

ubuntuPrettyName=$(cat /etc/os-release | grep PRETTY_NAME | cut -f2 -d'=')
ubuntuVersion=$(cat /etc/os-release | grep VERSION_ID | cut -f2 -d'=')
arch=$(uname -p)
echo "Detected $ubuntuPrettyName running on $arch."

isLatestOS=false
if [ $ubuntuVersion == "\"20.04\"" ] ; then
    isLatestOS=true
fi

isOldOS=false
if [ $ubuntuVersion \< "\"18.04\"" ] ; then
    isOldOS=true
fi

isX86=false
isArm=false
if [ $arch == "aarch64" ]; then
    isArm=true
fi
if [ $arch == "x86_64" ]; then
    isX86=true
fi
if [[ "$isX86" == false && "$isArm" == false ]]; then
    echo "Unsupported architecture. Only 'x86_64' and 'aarch64' are supported."
    exit 1
fi

if [[ "$isArm" == true && "$isLatestOS" == false ]]; then
    echo "Unsupported architecture & OS combination. Compiling on ARM requires Ubuntu 20.04 LTS or later."
    exit 1
fi

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


if [[ "$isLatestOS" == true ]]; then
    echo "Installing seqan v2.4.0 from Ubuntu repository."
    apt install -y libseqan2-dev
else
    # SeqAn 2.4.0 is not available in the older Ubuntu repositories.
    # Download it from GitHub, then install it.
    echo "Installing seqan 2.4.0 from Github."
    apt install -y curl
    curl -L https://github.com/seqan/seqan/releases/download/seqan-v2.4.0/seqan-library-2.4.0.deb \
        -o /dev/shm/seqan-library-2.4.0.deb
    apt install /dev/shm/seqan-library-2.4.0.deb
    rm /dev/shm/seqan-library-2.4.0.deb
fi


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

spoaBuildFlags="-Dspoa_generate_dispatch=ON"
if [[ "$isOldOS" == true || "$isArm" == true ]]; then
    spoaBuildFlags="-Dspoa_generate_dispatch=OFF -Dspoa_optimize_for_portability=OFF -Dspoa_optimize_for_native=OFF"
fi

# Build the shared library.
mkdir build
cd build
cmake ../spoa-v3.4.0 -DBUILD_SHARED_LIBS=ON $spoaBuildFlags
make -j all
make install

# Build the static library.
cd ..
mkdir build-static
cd build-static
cmake ../spoa-v3.4.0 -DBUILD_SHARED_LIBS=OFF $spoaBuildFlags
make -j all
make install
cd 
rm -rf $tmpDirectoryName


# Make sure the newly created libraries are immediately visible to the loader.
ldconfig


