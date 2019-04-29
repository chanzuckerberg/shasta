#!/usr/bin/python3

import shasta
import GetConfig
import argparse

# Get from the arguments the edge id (same as segment id).
parser = argparse.ArgumentParser(description='Write to a csv file assembly details for a single segment.')
parser.add_argument('edgeId', type=int)
edgeId = parser.parse_args().edgeId


# Read the config file.
config = GetConfig.getConfig()

# Create the Assembler.
a = shasta.Assembler()

# Set up the consensus caller.
a.setupConsensusCaller(config['Assembly']['consensusCaller'])

# Access what we need.
a.accessKmers()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.accessMarkerGraphVertexRepeatCounts()
a.accessMarkerGraphEdgeConsensus()
a.accessMarkerGraphCoverageData()
assembledSegment = a.assembleAssemblyGraphEdge(edgeId)

csv = open(str(edgeId) + '.csv', 'w')
for position in range(assembledSegment.size()):
    coverageData = assembledSegment.getCoverageData(position)
    csv.write('%i,' % position)
    csv.write('%s,' % assembledSegment.getBase(position))
    csv.write('%i,' % assembledSegment.getRepeatCount(position))
    for cd in coverageData:
       csv.write('%s%i%s %i,' % (cd.getBase(), cd.getRepeatCount(), cd.getStrand(), cd.getFrequency()))
    csv.write('\n')

