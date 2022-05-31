#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphEdges()
a.accessMode3AssemblyGraph()

path = a.createMode3PathGraph()


