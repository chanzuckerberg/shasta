#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphSequences()
a.writeGfa1('Assembly.gfa')



