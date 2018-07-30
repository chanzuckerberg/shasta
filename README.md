# chanzuckerberg/shasta

(See [here](https://chanzuckerberg.github.io/shasta/) for the complete documentation.)

It was [recently shown](https://www.nature.com/articles/nbt.4060) 
that de novo assembly for human genomes 
from [Oxford Nanopore](https://nanoporetech.com) reads is possible, 
but computationally very expensive and laborious. 
The aim of this project is to enable de novo assembly 
of human genomes on a routine/production basis 
and at reasonable computational cost. 
For this to be possible, de novo assembly must be:

* Fast: under one day elapsed time.
* Accurate: accuracy and other assembly metrics must be
comparable to those provided by existing tools.
Assembly accuracy should be limited by data quality,
not by assembly algorithms.
* Logistically simple and easy to run.

See [this presentation](https://chanzuckerberg.github.io/shasta/ShastaSlides-June2018-v2.pdf) 
for information on the computational approach currently being explored. 
In addition to Oxford Nanopore reads, these methods might also apply 
to long reads from other technologies, such as
the Pacific Biosciences DNA sequencing platforms.

This project is at an early stage. Its main output is currently
in the form of a *marker graph*, not assembled sequence. 
It does include functionality to extract and display a
local portion of the global marker graph, 
and to use it to assemble sequence using a simple-minded
algorithm.

Contributions of code, ideas, computational experiments, or documentation are welcome. 
To contribute, please use the standard GitHub Pull Request process. 
To facilitate and encourage contributions, the 
[guidelines for contributing](https://github.com/chanzuckerberg/shasta/blob/master/CONTRIBUTING.md)
are minimal.

Comments and criticism are also welcome. 
Please use the [Wiki](https://github.com/chanzuckerberg/shasta/wiki) 
or [Issues](https://github.com/chanzuckerberg/shasta/issues) 
sections of the repository as appropriate for these purposes.

For more detailed information,  see the complete documentation
[here](https://chanzuckerberg.github.io/shasta/).




#### Acknowledgment for file src/dset64.hpp

File `src/dset64.hpp` is a modified version of file `dset.h` from GitHub repository
[wjakob/dset](https://github.com/wjakob/dset) by Wenzel Jacob. 

See the LICENSE file for 
licensing information specific to dset64.hpp.

The code in `dset64.hpp` and `dset.h` implement the parallel 
algorithm described in 
*Wait-free Parallel Algorithms for the Union-Find Problem*
by Richard J. Anderson and Heather Woll,
STOC '91 Proceedings of the twenty-third annual ACM symposium on Theory of computing,
Pages 370-380.
It is used for efficient parallel
computation of the global marker graph.

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

