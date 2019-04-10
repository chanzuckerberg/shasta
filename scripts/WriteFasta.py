#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphSequences()
a.writeFasta('Assembly.fasta')



