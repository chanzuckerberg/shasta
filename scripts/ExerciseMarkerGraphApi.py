#!/usr/bin/python3

print(

"""
This scripts illustrates the use of the marker graph Python API.
It assumes that the run was already done (using RunAssembly.py)
and saved to disk (using SaveRun.py).
It must run from the run directory - that is, the directory
that contains dataOnDisk and DataOnDisk.

This uses vertex ids in the marker graph. To display these ids in the
browser, select "Detailed" and "Show vertex ids" when displaying
a local marker graph. The vertex ids displayed in this way
are the global vertex ids used by this script. 

""")


import shasta

a = shasta.Assembler(
    smallDataFileNamePrefix='dataOnDisk/',
    largeDataFileNamePrefix='DataOnDisk/')

a.accessReadsReadOnly()
a.accessMarkers()
a.accessMarkerGraphVertices()

while True:
    vertexIds = input('Enter a vertex id to find its children and parents\n' + 
        'or two vertex ids to find information about the edge that joins them:\n').split()
    if len(vertexIds)==1:
        vertexId = int(vertexIds[0]) 
        children = a.getGlobalMarkerGraphVertexChildren(
            globalMarkerGraphVertexId = vertexId)
        parents = a.getGlobalMarkerGraphVertexParents(
            globalMarkerGraphVertexId = vertexId)
        print('Start vertex:', vertexId)
        print('Parents:', parents)
        print('Children:', children)
    elif len(vertexIds)==2:
        parent = int(vertexIds[0])
        child = int(vertexIds[1])
        print('Parent:', parent)
        print('Child:', child)
        print('Not implemented.')

