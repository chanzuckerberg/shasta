#!/usr/bin/python3

import shasta

# Initialize the assembler and access what we need.
a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphReverseComplementEdge()
a.checkMarkerGraphIsStrandSymmetric()

