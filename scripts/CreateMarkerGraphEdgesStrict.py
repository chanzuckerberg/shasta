#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()


a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.createMarkerGraphEdgesStrict(
    minEdgeCoverage = int(config['MarkerGraph']['minEdgeCoverage']),
    minEdgeCoveragePerStrand = int(config['MarkerGraph']['minEdgeCoveragePerStrand']))


