#!/usr/bin/python3

import os
import psutil
import sys

"""

This copies to /dev/null all binary data needed by Mode2Assembly-B.py.
Invoke from the assembly directory, which contains the Data directory.

It can be used to get these data in cache, if enough memory is available.

If possible, for best performance, use 

sudo blockdev --setra 65536 /dev/xyz

then when done

sudo blockdev --setra 256 /dev/xyz

where xyz is the block device that contains the assembly directory.


"""

files = [
    'Info', 
    'Markers.toc', 'Markers.data', 
    'MarkerGraphVertexTable', 
    'MarkerGraphVertices.toc', 'MarkerGraphVertices.data',
    'GlobalMarkerGraphEdges',
    'GlobalMarkerGraphEdgeMarkerIntervals.toc', 'GlobalMarkerGraphEdgeMarkerIntervals.data',
    'GlobalMarkerGraphEdgesBySource.toc', 'GlobalMarkerGraphEdgesBySource.data',
    'GlobalMarkerGraphEdgesByTarget.toc', 'GlobalMarkerGraphEdgesByTarget.data',
    'MarkerGraphEdgesConsensus.toc', 'MarkerGraphEdgesConsensus.data',
    'MarkerGraphEdgesConsensusOverlappingBaseCount']
    
gb = 1024 * 1024 * 1024
totalSize = 0
for file in files:
    totalSize = totalSize + os.path.getsize("Data/%s" % file)
print('Total size of the required binary data:', int(round(totalSize / gb)), 'GB')   

memory =  psutil.virtual_memory().total
print('Total physical memory:', int(round(memory / gb)), 'GB')


if totalSize > memory:
    sys.exit('Not enough memory.')
    

for file in files:
    os.system('cp Data/%s /dev/null' % file)
    
print('Done')





