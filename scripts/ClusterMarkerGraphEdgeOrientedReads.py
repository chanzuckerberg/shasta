#!/usr/bin/python3

import shasta
import GetConfig
import argparse

parser = argparse.ArgumentParser(
    description='Cluster the oriented reads of a marker graph edge based on theri sequence.')
parser.add_argument('edgeId', type=int)
arguments = parser.parse_args()
edgeId = arguments.edgeId



config = GetConfig.getConfig()

a = shasta.Assembler()

a.accessMarkers()
a.accessMarkerGraphEdges(True, True)

a.clusterMarkerGraphEdgeOrientedReads(edgeId)


