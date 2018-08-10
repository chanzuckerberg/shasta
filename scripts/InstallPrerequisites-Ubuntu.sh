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

