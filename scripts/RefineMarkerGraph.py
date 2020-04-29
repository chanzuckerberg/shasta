#!/usr/bin/python3

import shasta

# This script destroys the original marker graph.
# To recreate it and bring it back to its state
# after transitive reduction, use 
# CreateMarkerGraphAndTransitiveReduction.py

a = shasta.Assembler()
a.accessMarkerGraphVertices(readWriteAccess = True)
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphReverseComplementEdge()
a.refineMarkerGraph()

