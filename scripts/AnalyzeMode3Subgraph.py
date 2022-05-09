#!/usr/bin/python3

import shasta

segmentIds = [int(token) for token in input('Enter segment ids on one line: ').split()] 

a = shasta.Assembler()
a.accessMode3AssemblyGraph()
a.analyzeMode3Subgraph(segmentIds)


