#!/usr/bin/python3

import shasta
import GetConfig
import ast



# Create the Assembler.
a = shasta.Assembler()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.findAssemblyGraphBubbles()



