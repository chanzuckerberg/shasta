#!/bin/bash

# This script, which must be run with root privilege,
# installs the version of the CUDA toolkit
# that is required to run Shasta with GPU acceleration
# on Ubuntu 18.04.
# This replicates the instructions given here:
# https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&target_distro=Ubuntu&target_version=1804&target_type=deblocal

# Create a temporary directory.
tmpDirectoryName=$(mktemp --directory --tmpdir)
echo $tmpDirectoryName
cd $tmpDirectoryName

apt-get update
apt install -y wget
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"
apt-get update
apt-get -y install cuda

# Remove the temporary directory.
cd 
rm -rf $tmpDirectoryName

