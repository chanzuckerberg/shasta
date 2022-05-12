#!/usr/bin/python3

import shasta

(segmentId, direction) = (int(token) for token in input('Enter segment id and direction on one line: ').split()) 

a = shasta.Assembler()
a.accessMode3AssemblyGraph()

path = []
a.createMode3AssemblyPath(segmentId, direction, path)

print('Found the following path: ', *path)
