This directory contains experimental code used to run portions of Shasta on a GPU. 

This code is not built by default. To build it, use the following when running cmake:

-DBUILD_FOR_GPU=ON

All GPU code that is not in this directory should be under 

#ifdef SHASTA_BUILD_FOR_GPU


