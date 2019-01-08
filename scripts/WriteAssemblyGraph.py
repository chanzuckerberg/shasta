#!/usr/bin/python3

import shasta
import sys

if not len(sys.argv) == 2:
     raise Exception('Call with one argument, the file name.')
fileName = sys.argv[1];

a = shasta.Assembler()
a.accessAssemblyGraphEdges()
a.writeAssemblyGraph(fileName)



