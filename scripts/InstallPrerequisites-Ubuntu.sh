#!/bin/bash

# This installs all packages needed to build Shasta or to run
# the dynamic executable. The static executable has no prerequisites.

minimalInstall=false
for arg in "$@"
do
    case $arg in
        --minimal)
        minimalInstall=true
        shift
        ;;
        *)
        echo "Usage: /path/to/InstallPrerequisites-Ubuntu.sh [--minimal]"
        exit 1
        ;;
    esac
done

ubuntuPrettyName=$(cat /etc/os-release | grep PRETTY_NAME | cut -f2 -d'=')
ubuntuVersion=$(cat /etc/os-release | grep VERSION_ID | cut -f2 -d'=')
arch=$(uname -p)
echo "Detected $ubuntuPrettyName running on $arch."

isLatestOS=true
if [ $ubuntuVersion \< "\"20.04\"" ] ; then
    isLatestOS=false
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

if [ "$isLatestOS" == false ]; then
    echo "Unsupported OS. Building Shasta requires Ubuntu 20.04 LTS or later."
    exit 1
fi

apt-get update
apt-get install --yes \
    git \
    g++ \
    make \
    cmake \
    curl \
    libboost-system-dev \
    libboost-program-options-dev \
    libboost-chrono-dev \
    libpng-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    libseqan2-dev

if [ "$minimalInstall" == false ]; then
    # Install packages required for the HTTP server and Python-C++ bindings.
    apt-get install --yes \
        libboost-graph-dev \
        ncbi-blast+ \
        graphviz \
        gnuplot \
        python3-dev \
        python3-pip
        pip3 install pybind11==2.5.0
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
if [[ "$isArm" == true ]]; then
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

echo "=============================="
echo " Relevant software versions"
echo "=============================="
if [ "$minimalInstall" == false ]; then
    echo "$(apt list git g++ cmake curl libseqan2-dev libboost-system-dev libpng-dev graphviz gnuplot ncbi-blast+ python3 python3-pip 2>/dev/null)"
else
    echo "$(apt list git g++ cmake curl libseqan2-dev libboost-system-dev 2>/dev/null)"
fi
echo "=============================="

# Make sure the newly created libraries are immediately visible to the loader.
ldconfig
