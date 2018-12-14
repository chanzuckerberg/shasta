#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessKmers()
a.accessReadsReadOnly()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphConnectivity()
a.accessAssemblyGraphVertices()
a.setupConsensusCaller(config['Assembly']['consensusCaller'])
a.assemble()



