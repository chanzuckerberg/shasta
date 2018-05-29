#!/usr/bin/python3

import os
import sys



if len(sys.argv) != 4:
    print('Creates a local subgraph around a start vertex.\n'
          'Invoke with three arguments:\n',
          '- Name of the input dot file.\n',
          '- Name of the start vertex.\n',
          '- Maximum distance.')
          
inputFileName = sys.argv[1]          
vertexName = sys.argv[2]          
maxDistance = sys.argv[3]  

# Create an intermediate file that contains distances from the start vertex
# (vertex attribute "dist").
intermediateFileName =  'WithDistances-' + inputFileName
os.system(' '.join(['dijkstra -a', '"' + vertexName + '"', inputFileName, '>', intermediateFileName]))

# Only keep vertices that are close enough.
outputFileName = 'LocalSubgraph-' + inputFileName
os.system('gvpr -i \'N[dist<%s.]\' %s -o %s' % (maxDistance, intermediateFileName, outputFileName))

os.remove(intermediateFileName)

