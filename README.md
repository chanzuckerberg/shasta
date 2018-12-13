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



#### Key computational ideas

See [this somewhat obsolete presentation](https://chanzuckerberg.github.io/shasta/ShastaSlides-June2018-v2.pdf) 
for information on the computational approach currently being explored. 
Key computational ideas under considerations include:

* Use a fixed subset of all k-mers (with short k, typically 10), called markers,
as a concise representation of reads useful to find pairs of overlapping reads 
and to efficiently compute alignments.
* Use a [MinHash](https://en.wikipedia.org/wiki/MinHash) approach
on this concise representation to find pairs
of potentially overlapping reads.
* To check if two reads actually overlap, compute a marker alignment, rather than a 
complete alignment.
* Use a [run-length encoding](https://en.wikipedia.org/wiki/Run-length_encoding) of reads, 
as first proposed 
[here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0046679) for
error correction of Pacific Biosciences reads.
The run-length representation of reads in insensitive to one of the most
common types of errors in Nanopore reads, namely errors in the length of a homopolymer run.

In addition to Oxford Nanopore reads, these methods might also apply 
to long reads from other technologies, such as
the Pacific Biosciences DNA sequencing platforms.



#### Current status

This project is at an early stage. Its main output is currently
in the form of a *marker graph*, not assembled sequence. 
It does include functionality to extract and display a
local portion of the global marker graph, 
and to use it interactively for local assembly.

This limited computation currently runs in a few hours for
a typical human genome with coverage 40x, using a single
machine with sufficient memory and 32 cores.



#### Acknowledgments

The Shasta software uses various external software packages.
See [https://chanzuckerberg.github.io/shasta/Acknowledgments.html](here) for more information.



#### Contributions and comments

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






