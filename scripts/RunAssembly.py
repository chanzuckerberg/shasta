#!/usr/bin/python3

import shasta
import GetConfig
import ast
import os
import sys


def verifyConfigFiles(parentDirectory=""):
    # Generate absolute paths to the files that will be created
    dataPath = os.path.abspath(os.path.join(parentDirectory, "Data"))
    confPath = os.path.abspath(os.path.join(parentDirectory, "shasta.conf"))

    # Check that we have what we need.
    if not os.path.lexists(dataPath):
        raise Exception('Missing: Data. Use SetupRunDirectory.py to set up the run directory.')
    if not os.path.lexists(confPath):
        raise Exception('Missing: configuration file shasta.conf. Sample available in shasta-install/conf.')


def verifyFastaFiles(fastaFileNames):
    # Get from the arguments the list of input fasta files and check that they all exist.
    helpMessage = "Invoke passing as arguments the names of the input Fasta files."
    if len(sys.argv)==1:
        print(helpMessage)
        exit(1)

    for fileName in fastaFileNames:
        if not os.path.lexists(fileName):
            raise Exception('Input file %s not found.' % fileName)


def initializeAssembler(config, fastaFileNames):
    # Create the Assembler.
    a = shasta.Assembler(createNew=True)    
    return a


def runAssembly(a, config, fastaFileNames):
    # Set up the consensus caller.
    a.setupConsensusCaller(config['Assembly']['consensusCaller'])
    
    # Figure out if we should use marginPhase, and if so set it up.
    useMarginPhase = ast.literal_eval(config['Assembly']['useMarginPhase'])
    if useMarginPhase:
        a.setupMarginPhase()
    
    # Read the input fasta files.
    for fileName in fastaFileNames:  
        print('Reading input file', fileName, flush=True) 
        a.addReadsFromFasta(
            fileName = fileName, 
            minReadLength = int(config['Reads']['minReadLength']))

    # Initialize read flags.
    a.initializeReadFlags()            
    
    # Create a histogram of read lengths.
    a.histogramReadLength(fileName="ReadLengthHistogram.csv")
    
    # Randomly select the k-mers that will be used as markers.
    a.randomlySelectKmers(
        k = int(config['Kmers']['k']), 
        probability = float(config['Kmers']['probability']))
        
    # Find the markers in the reads.
    a.findMarkers()
    
    # Flag palindromic reads.
    # These wil be excluded from further processing.
    a.flagPalindromicReads(
        maxSkip = int(config['Reads']['palindromicReads.maxSkip']),
        maxMarkerFrequency = int(config['Reads']['palindromicReads.maxMarkerFrequency']),
        alignedFractionThreshold = float(config['Reads']['palindromicReads.alignedFractionThreshold']),
        nearDiagonalFractionThreshold = float(config['Reads']['palindromicReads.nearDiagonalFractionThreshold']),
        deltaThreshold = int(config['Reads']['palindromicReads.deltaThreshold']))
        
    # Find alignment candidates.
    a.findAlignmentCandidatesLowHash(
        m = int(config['MinHash']['m']), 
        hashFraction = float(config['MinHash']['hashFraction']),
        minHashIterationCount = int(config['MinHash']['minHashIterationCount']), 
        maxBucketSize = int(config['MinHash']['maxBucketSize']),
        minFrequency = int(config['MinHash']['minFrequency']))
    """
    # Old MinHash code to find alignment candidates. 
    # If using this, make sure to set MinHash.minHashIterationCount
    # to a suitable value (100 was the defaulf when this code was active).
    a.findAlignmentCandidatesMinHash(
        m = int(config['MinHash']['m']), 
        minHashIterationCount = int(config['MinHash']['minHashIterationCount']), 
        maxBucketSize = int(config['MinHash']['maxBucketSize']),
        minFrequency = int(config['MinHash']['minFrequency']))
    """
    
    # Compute alignments.
    a.computeAlignments(
        maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
        maxSkip = int(config['Align']['maxSkip']),
        minAlignedMarkerCount = int(config['Align']['minAlignedMarkerCount']),
        maxTrim = int(config['Align']['maxTrim']))
        
    # Create the read graph.
    a.createReadGraph(
        maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']),
        maxTrim = int(config['Align']['maxTrim']))

    # Flag read graph edges that cross strands.
    a.flagCrossStrandReadGraphEdges()
    
    # Flag chimeric reads.
    a.flagChimericReads(
        maxChimericReadDistance = int(config['ReadGraph']['maxChimericReadDistance']))
    a.computeReadGraphConnectedComponents(
        minComponentSize = int(config['ReadGraph']['minComponentSize']))
    
    # Create vertices of the marker graph.
    a.createMarkerGraphVertices(
        maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
        maxSkip = int(config['Align']['maxSkip']),
        minCoverage = int(config['MarkerGraph']['minCoverage']),
        maxCoverage = int(config['MarkerGraph']['maxCoverage']))
    a.findMarkerGraphReverseComplementVertices()
    
    # Create edges of the marker graph.
    a.createMarkerGraphEdges()
    a.findMarkerGraphReverseComplementEdges()
    
    # Approximate transitive reduction.
    a.transitiveReduction(
        lowCoverageThreshold = int(config['MarkerGraph']['lowCoverageThreshold']),
        highCoverageThreshold = int(config['MarkerGraph']['highCoverageThreshold']),
        maxDistance = int(config['MarkerGraph']['maxDistance']),
        edgeMarkerSkipThreshold = int(config['MarkerGraph']['edgeMarkerSkipThreshold'])
        )
    
    # Prune the strong subgraph of the marker graph.
    a.pruneMarkerGraphStrongSubgraph(
        iterationCount = int(config['MarkerGraph']['pruneIterationCount']))
    
    # Simplify the marker graph to remove bubbles and superbubbles.
    # The maxLength parameter controls the maximum number of markers
    # for a branch to be collapsed during each iteration.
    a.simplifyMarkerGraph(
        maxLength = [int(s) for s in config['MarkerGraph']['simplifyMaxLength'].split(',')],
        debug = False)
    
    # Create the assembly graph.
    a.createAssemblyGraphEdges()
    a.createAssemblyGraphVertices()
    a.writeAssemblyGraph("AssemblyGraph-Final.dot")
    
    # Compute optimal repeat counts for each vertex of the marker graph.
    a.assembleMarkerGraphVertices()
    
    # Optionally compute coverage data for marker graph vertices.
    storeCoverageData = ast.literal_eval(config['Assembly']['storeCoverageData'])
    if storeCoverageData:
        a.computeMarkerGraphVerticesCoverageData()
    
    # Compute consensus sequence for marker graph edges to be used for assembly.
    a.assembleMarkerGraphEdges(
        markerGraphEdgeLengthThresholdForConsensus =
        int(config['Assembly']['markerGraphEdgeLengthThresholdForConsensus']),
        useMarginPhase = useMarginPhase,
        storeCoverageData = storeCoverageData)
    
    # Use the assembly graph for global assembly.
    a.assemble()
    
    a.computeAssemblyStatistics()
    a.writeGfa1('Assembly.gfa')
    a.writeFasta('Assembly.fasta')


def main():
    # Parse arguments
    fastaFileNames = sys.argv[1:]

    # Ensure prerequisite files are present
    verifyConfigFiles()
    verifyFastaFiles(fastaFileNames=fastaFileNames)

    # Read the config file.
    config = GetConfig.getConfig()

    # Initialize Assembler object
    assembler = initializeAssembler(config=config, fastaFileNames=fastaFileNames)

    # Run with user specified configuration and input files
    runAssembly(config=config, fastaFileNames=fastaFileNames, a=assembler)


if __name__ == "__main__":
    main()
