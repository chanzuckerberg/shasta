#!/usr/bin/python3

import shasta

segmentIds = [200, 300, 400]

a = shasta.Assembler()
a.accessMode3AssemblyGraph()
a.analyzeMode3Subgraph(segmentIds)


