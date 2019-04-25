#!/usr/bin/python3

print(

"""
This scripts illustrates the use of Python-callable function 
Assembler.getLocalAssemblyPath.
It must run from the run directory - that is, the directory
that contains dataOnDisk and DataOnDisk.

This uses vertex ids in the marker graph. To display these ids in the
browser, select "Detailed" and "Show vertex ids" when displaying
a local marker graph. The vertex ids displayed in this way
are the global vertex ids used by this script. 

""")


import shasta

# Create the Assembler object and access what we need.
a = shasta.Assembler(
    largeDataFileNamePrefix='DataOnDisk/')
a.accessMarkers()
a.accessMarkerGraphVertices()

# Interactive loop.
while True:

    # Get the start vertex.
    tokens = input('Specify the start vertex in one of two ways:\n' + 
        '- By its vertex id.\n'
        '- By its read id, strand, marker ordinal (separated by spaces).\n').split()
        
    # If we have one token, interpret is as the vertex id of the start vertex.
    if len(tokens) == 1:
        startVertexId = int(tokens[0])
        
    # If we have three tokens, interpret them as read id, strand, marker ordinal.
    elif len(tokens) == 3:
        readId = int(tokens[0])
        strand = int(tokens[1])
        ordinal = int(tokens[2])
        startVertexId = a.getGlobalMarkerGraphVertex(readId, strand, ordinal)
        if startVertexId == shasta.invalidCompressedGlobalMarkerGraphVertexId:
            print('There is no vertex associated with this marker.')
            continue
    else:
        print('Nothing done, try again.')
        continue;
    print('Start vertex id is ', startVertexId)
        
    # Get the maximum distance.
    maxDistance = int(input('Enter the maximum distance from the start vertex:\n'))
    print('Maximum distance from start vertex is', maxDistance)
    
    # Create the local marker graph and find its local assembly path.
    path = a.getLocalAssemblyPath(startVertexId=startVertexId, maxDistance=maxDistance)
    
    print('Assembly path:', path)

