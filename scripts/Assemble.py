#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphConnectivity()
a.accessAssemblyGraphVertices()
a.setupConsensusCaller('SimpleConsensusCaller')
a.assemble()



