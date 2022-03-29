#!/usr/bin/python3

import shasta
import GetConfig
import argparse

parser = argparse.ArgumentParser(
    description='Cluster the oriented reads of a marker graph edge based on their sequence.')
parser.add_argument('edgeId', type=int)
parser.add_argument('errorRateThreshold', type=float)
arguments = parser.parse_args()
edgeId = arguments.edgeId
errorRateThreshold = arguments.errorRateThreshold



config = GetConfig.getConfig()

a = shasta.Assembler()

a.accessMarkers()
a.accessMarkerGraphEdges(True, True)

a.clusterMarkerGraphEdgeOrientedReads(
    edgeId = edgeId, 
    errorRateThreshold = errorRateThreshold,
    debug = True)


