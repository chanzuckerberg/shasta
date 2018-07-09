# chanzuckerberg/shasta

It was [recently shown](https://www.nature.com/articles/nbt.4060) 
that de novo assembly for human genomes 
from [Oxford Nanopore](https://nanoporetech.com) reads is possible, 
but computationally expensive and laborious. 
The aim of this project is to make de novo assembly 
of human genomes possible on a routine/production basis 
and at reasonable computational cost. 
For this to be possible, de novo assembly must be:

* Fast: under one day elapsed time.
* Accurate: accuracy and other assembly metrics comparable to those provided by existing tools.
* Logistically simple and easy to run.

This project is at an early stage, and does not yet 
produce assembled sequence as its final output. 
Contributions of code, ideas, computational experiments, or documentation are welcome. 
To contribute, please use the standard GitHub Pull Request process. 
To facilitate and encourage contributions, the guidelines for contributing are minimal.

Comments and criticism are also welcome. 
Please use the [Wiki](https://github.com/chanzuckerberg/shasta/wiki) 
or [Issues](https://github.com/chanzuckerberg/shasta/issues) 
sections of the repository as appropriate for these purposes.

See [this presentation](docs/ShastaSlides-June2018-v2.pdf) 
for information on the computational approach selected. 
In addition to Oxford Nanopore reads, these methods might also apply 
to other long reads with high error rates 
such as those created by the Pacific Biosciences DNA sequencing platforms.



#### Acknowledgment for file src/dset64.hpp

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
which brings the maximum number of items to 2<sup>64</sup>.

Many thanks to Wenzel Jacob for making his implementation
available as open source software, and for providing 
helpful information during the conversion process to 64 bits.

