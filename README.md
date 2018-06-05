# Nanopore2

De novo assembly from Nanopore sequencing data.
This repository is at an early stage and
no additional documentation is available.


### Acknowledgment

File src/dset64.hpp is a modified version of file dset.h from GitHub repository
wjakob/dset by Wenzel Jacob. See the LICENSE file for 
licensing information specific to dset64.hpp.

Both files implement the parallel 
algorithm described in 
*Wait-free Parallel Algorithms for the Union-Find Problem*
by Richard J. Anderson and Heather Woll,
STOC '91 Proceedings of the twenty-third annual ACM symposium on Theory of computing,
Pages 370-380.

See [this Wikipedia article](https://en.wikipedia.org/wiki/Disjoint-set_data_structure)
for more information on the sequential version of the algorithm.

The original implementation by Wenzel Jacob uses 64-bit
atomic primitives, and implements the union-find 
algorithm for 32-bit item ids, which allows up to 
2<sup>32</sup> items. The modified version in src/dset64.hpp
uses 128-bit primitives for 64-bit item ids,
which brings the maximum number of items to 2^64^.

Many thanks to Wenzel Jacob for making his implementation
available as open source software, and for providing 
helpful information during the conversion process to 64 bits.

