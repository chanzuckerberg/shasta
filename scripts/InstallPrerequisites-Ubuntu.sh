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
    libboost-graph-dev \
    libboost-chrono-dev \
    libpng-dev \
    libblas-dev \
    liblapack-dev \
    gfortran \
    libseqan2-dev

if [ "$minimalInstall" == false ]; then
    # Install packages required for the HTTP server and Python-C++ bindings.
    apt-get install --yes \
        ncbi-blast+ \
        graphviz \
        gnuplot \
        python3-dev \
        python3-pip
        pip3 install pybind11==2.8.1
fi



# The spoa library is available in the stable Ubuntu repository, but
# without the static version.
# So we have to build it from source.

# Create a temporary directory.
tmpDirectoryName=$(mktemp --directory --tmpdir)
echo $tmpDirectoryName
cd $tmpDirectoryName

# Get the code.
curl -L https://github.com/rvaser/spoa/archive/refs/tags/4.0.8.tar.gz \
    -o 4.0.8.tar.gz
tar -xvf 4.0.8.tar.gz

# The spoa dispatcher feature selects code at run time based on available hardware features,
# which can improve performance.
# However, in spoa v4.0.8 it introduces two additional dependencies:
# - USCiLab/cereal
# - google/cpu_features
# To avoid these additional dependencies, we turn off the dispatcher feature for now.
# We could turn it back on if we see significant performance degradation in this area.
spoaBuildFlags="-Dspoa_generate_dispatch=ON"
if [[ "$isArm" == true ]]; then
    spoaBuildFlags="-Dspoa_generate_dispatch=OFF -Dspoa_optimize_for_portability=OFF -Dspoa_optimize_for_native=OFF"
fi
# Per the above comment, turn off the dispatcher feature for now.
spoaBuildFlags="-DCMAKE_BUILD_TYPE=Release"



# Build the shared library.
mkdir build
cd build
cmake ../spoa-4.0.8 -DBUILD_SHARED_LIBS=ON $spoaBuildFlags
make -j all
make install

# Build the static library.
cd ..
mkdir build-static
cd build-static
cmake ../spoa-4.0.8 -DBUILD_SHARED_LIBS=OFF $spoaBuildFlags
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
