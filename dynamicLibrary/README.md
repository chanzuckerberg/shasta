# Directory shasta /dynamicLibrary

This directory builds the Shasta dynamic library `shasta.so`.
This library is used in three ways:

* It is linked in by the shasta dynamic executable `shastaDynamic`.
* It can be imported by a python script via `import shasta` to provide Shasta Python bindings.
* It can be statically linked in by other C++ code outside Shasta that uses Shasta as a library.
