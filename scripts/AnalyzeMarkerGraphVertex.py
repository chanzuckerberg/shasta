#!/usr/bin/python3

import argparse
import shasta

# Get the id of the vertex we want to analyze.
parser = argparse.ArgumentParser(description=
    'Analyze a vertex of the marker graph.')   
parser.add_argument('vertexId', type=int)
arguments = parser.parse_args()
vertexId = arguments.vertexId

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessDirectedReadGraphReadOnly()
a.accessConflictReadGraph()
a.analyzeMarkerGraphVertex(vertexId)
