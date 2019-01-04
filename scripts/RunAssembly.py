#!/usr/bin/python3

import shasta
import GetConfig
import os
import sys

# Check that we have what we need.
if not os.path.lexists('Data'):
    raise Exception('Missing: Data. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('threadLogs'):
    raise Exception('Missing: threadLogs. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('shasta.conf'):
    raise Exception('Missing: configuration file shasta.conf. Sample available in shasta-install/conf.')

# Get from the arguments the list of input fasta files and check that they all exist.
helpMessage = "Invoke passing as arguments the names of the input Fasta files."
if len(sys.argv)==1:
    print(helpMessage)
    exit(1)
fastaFileNames = sys.argv [1:];
for fileName in fastaFileNames:  
    if not os.path.lexists(fileName):
        raise Exception('Input file %s not found.' % fileName)

# Read the config file.
config = GetConfig.getConfig()

# Create the assembler.
useRunLengthReadsString = config['Reads']['useRunLengthReads']
if useRunLengthReadsString == 'True':
    useRunLengthReads = True
elif useRunLengthReadsString == 'False':
    useRunLengthReads = False
else:
    raise RuntimeError("Configuration parameter useRunLengthReads in section Reads must be True or False.")
a = shasta.Assembler(useRunLengthReads = useRunLengthReads)

# Read the input fasta files.
a.accessReadsReadWrite();
a.accessReadNamesReadWrite();
for fileName in fastaFileNames:  
    print('Reading input file', fileName, flush=True) 
    a.addReadsFromFasta(
        fileName = fileName, 
        minReadLength = int(config['Reads']['minReadLength']))

# Create a histogram of read lengths.
a.histogramReadLength(fileName="ReadLengthHistogram.csv")

# Randomly select the k-mers that will be used as markers.
a.randomlySelectKmers(
    k = int(config['Kmers']['k']), 
    probability = float(config['Kmers']['probability']))
    
# Find the markers in the reads.
a.findMarkers()

# Run MinHash to find pairs of reads that may overlap.
a.findAlignmentCandidates(
    m = int(config['MinHash']['m']), 
    minHashIterationCount = int(config['MinHash']['minHashIterationCount']), 
    log2MinHashBucketCount = int(config['MinHash']['log2MinHashBucketCount']),
    maxBucketSize = int(config['MinHash']['maxBucketSize']),
    minFrequency = int(config['MinHash']['minFrequency']))

# Compute alignments.
a.computeAlignments(
    maxVertexCountPerKmer = int(config['Align']['maxVertexCountPerKmer']),
    maxSkip = int(config['Align']['maxSkip']),
    minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
    maxTrim = int(config['Align']['maxTrim']))
    
# Create the read graph.
a.createReadGraph(maxTrim = int(config['Align']['maxTrim']))

# Flag chimeric reads.
a.flagChimericReads(
    maxChimericReadDistance = int(config['ReadGraph']['maxChimericReadDistance']))

# Create vertices of the marker graph.
a.createMarkerGraphVertices(
    maxVertexCountPerKmer = int(config['Align']['maxVertexCountPerKmer']),
    maxSkip = int(config['Align']['maxSkip']),
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']))

# Create edges of the marker graph.
a.createMarkerGraphEdges()

a.flagMarkerGraphWeakEdges(
    lowCoverageThreshold = int(config['MarkerGraph']['lowCoverageThreshold']),
    highCoverageThreshold = int(config['MarkerGraph']['highCoverageThreshold']),
    maxDistance = int(config['MarkerGraph']['maxDistance']),
    )

# Prune the strong subgraph of the marker graph.
a.pruneMarkerGraphStrongSubgraph(
    iterationCount = int(config['MarkerGraph']['pruneIterationCount']))

# Create the assembly graph.
a.createAssemblyGraphEdges()
a.createAssemblyGraphVertices()

# Use the assembly graph for global assembly.
a.setupConsensusCaller(config['Assembly']['consensusCaller'])
a.assemble()
a.computeAssemblyStatistics()
a.writeGfa1('Assembly.gfa')




