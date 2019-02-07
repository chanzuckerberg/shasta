#!/usr/bin/python3

from SetupSmallRunDirectory import setupSmallRunDirectory, verifyDirectoryFiles
from RunAssembly import verifyConfigFiles, verifyFastaFiles, runAssembly, initializeAssembler
import configparser

from datetime import datetime
from shutil import copyfile
import subprocess
import argparse
import sys
import os


def getDatetimeString():
    """
    Generate a datetime string. Useful for making output folders names that never conflict.
    """
    now = datetime.now()
    now = [now.year, now.month, now.day, now.hour, now.minute, now.second, now.microsecond]
    datetimeString = "_".join(list(map(str, now)))

    return datetimeString


def ensureDirectoryExists(directoryPath, i=0):
    """
    Recursively test directories in a directory path and generate missing directories as needed
    :param directoryPath:
    :return:
    """
    if i > 3:
        print("WARNING: generating subdirectories of depth %d, please verify path is correct: %s" % (i, directoryPath))

    if not os.path.exists(directoryPath):
        try:
            os.mkdir(directoryPath)

        except FileNotFoundError:
            ensureDirectoryExists(os.path.dirname(directoryPath), i=i + 1)

            if not os.path.exists(directoryPath):
                os.mkdir(directoryPath)
                

def overrideDefaultConfig(config, args):
    """
    Check all the possible params to see if the user provided an override value, and add any overrides
    to their appropriate location in the config dictionary
    """
    if args.minReadLength is not None:
        config["Reads"]["minReadLength"] = str(args.minReadLength)
        
    if args.k is not None:
        config["Kmers"]["k"] = str(args.k)
        
    if args.probability is not None:
        config["Kmers"]["probability"] = str(args.probability)
        
    if args.m is not None:
        config["MinHash"]["m"] = str(args.m)
        
    if args.minHashIterationCount is not None:
        config["MinHash"]["minHashIterationCount"] = str(args.minHashIterationCount)
        
    if args.maxBucketSize is not None:
        config["MinHash"]["maxBucketSize"] = str(args.maxBucketSize)
        
    if args.minFrequency is not None:
        config["MinHash"]["minFrequency"] = str(args.minFrequency)
        
    if args.maxSkip is not None:
        config["Align"]["maxSkip"] = str(args.maxSkip)
        
    if args.maxMarkerFrequency is not None:
        config["Align"]["maxMarkerFrequency"] = str(args.maxMarkerFrequency)
        
    if args.minAlignedMarkerCount is not None:
        config["Align"]["minAlignedMarkerCount"] = str(args.minAlignedMarkerCount)
        
    if args.maxTrim is not None:
        config["Align"]["maxTrim"] = str(args.maxTrim)
        
    if args.minComponentSize is not None:
        config["ReadGraph"]["minComponentSize"] = str(args.minComponentSize)
        
    if args.maxChimericReadDistance is not None:
        config["ReadGraph"]["maxChimericReadDistance"] = str(args.maxChimericReadDistance)
        
    if args.minCoverage is not None:
        config["MarkerGraph"]["minCoverage"] = str(args.minCoverage)
        
    if args.maxCoverage is not None:
        config["MarkerGraph"]["maxCoverage"] = str(args.maxCoverage)
        
    if args.lowCoverageThreshold is not None:
        config["MarkerGraph"]["lowCoverageThreshold"] = str(args.lowCoverageThreshold)
        
    if args.highCoverageThreshold is not None:
        config["MarkerGraph"]["highCoverageThreshold"] = str(args.highCoverageThreshold)
        
    if args.maxDistance is not None:
        config["MarkerGraph"]["maxDistance"] = str(args.maxDistance)
        
    if args.pruneIterationCount is not None:
        config["MarkerGraph"]["pruneIterationCount"] = str(args.pruneIterationCount)
        
    if args.shortCycleLengthThreshold is not None:
        config["MarkerGraph"]["shortCycleLengthThreshold"] = str(args.shortCycleLengthThreshold)
        
    if args.bubbleLengthThreshold is not None:
        config["MarkerGraph"]["bubbleLengthThreshold"] = str(args.bubbleLengthThreshold)
        
    if args.superBubbleLengthThreshold is not None:
        config["MarkerGraph"]["superBubbleLengthThreshold"] = str(args.superBubbleLengthThreshold)
        
    if args.markerGraphEdgeLengthThresholdForConsensus is not None:
        config["Assembly"]["markerGraphEdgeLengthThresholdForConsensus"] = str(args.markerGraphEdgeLengthThresholdForConsensus)
        
    if args.consensusCaller is not None:
        config["Assembly"]["consensusCaller"] = str(args.consensusCaller) + "ConsensusCaller"
        
    if args.useMarginPhase is not None:
        config["Assembly"]["useMarginPhase"] = str(args.useMarginPhase)

    return config


def main(readsSequencePath, outputParentDirectory, args):    
    if not os.path.exists(readsSequencePath):
        raise Exception("ERROR: input file not found: %s" % readsSequencePath)

    # Make sure given sequence file path is absolute, because CWD will be changed later
    readsSequencePath = os.path.abspath(readsSequencePath)

    # Generate output directory to run shasta in
    outputDirectoryName = "run_" + getDatetimeString()
    outputDirectory = os.path.join(outputParentDirectory, outputDirectoryName)
    ensureDirectoryExists(outputDirectory)
    
    # Setup run directory according to SetupSmallRunDirectory.py
    verifyDirectoryFiles(parentDirectory=outputDirectory)
    setupSmallRunDirectory(outputDirectory)

    # Locate path of default configuration files relative to this script's "binary" file
    scriptPath = os.path.abspath(os.path.dirname(__file__))
    confDirectory = os.path.join(os.path.dirname(scriptPath), "conf")

    defaultConfFilename = "shasta.conf"
    defaultConfPath = os.path.join(confDirectory, defaultConfFilename)
    localConfPath = os.path.join(outputDirectory, "shasta.conf")

    # Parse config file to fill in default parameters
    config = configparser.ConfigParser()
    if not config.read(defaultConfPath):
        raise Exception("Error reading config file %s." % defaultConfPath)

    # Check if any params were specified by user and override the default config
    config = overrideDefaultConfig(config, args)

    # Write updated config file to output directory so RunAssembly.py can be called as a separate process
    with open(localConfPath, "w") as file:
        config.write(file)

    # Add bayesian params file to the output directory if needed
    if args.consensusCaller == "SimpleBayesian":
        defaultMatrixPath = os.path.join(confDirectory, "SimpleBayesianConsensusCaller-1.csv")
        localMatrixPath = os.path.join(outputDirectory, "SimpleBayesianConsensusCaller.csv")
        copyfile(defaultMatrixPath, localMatrixPath)

    # Add marginphase params file to the output directory if needed
    if args.useMarginPhase:
        defaultParamsPath = os.path.join(confDirectory, "MarginPhase-allParams.np.json")
        localParamsPath = os.path.join(outputDirectory, "MarginPhase.json")
        copyfile(defaultParamsPath, localParamsPath)

    # Ensure prerequisite files are present
    verifyConfigFiles(parentDirectory=outputDirectory)
    verifyFastaFiles(fastaFileNames=[readsSequencePath])

    # Set current working directory to the output dir
    os.chdir(outputDirectory)
    
    executablePath = os.path.join(scriptPath, "RunAssembly.py")
    arguments = [executablePath, readsSequencePath]

    # processHandler.launchProcess(Data=Data, largePagesMountPoint=largePagesMountPoint, arguments=arguments)
    assemblyProcess = subprocess.Popen(arguments)


def stringAsBool(s):
    s = s.lower()
    boolean = None

    if s in {"t", "true", "1", "y", "yes"}:
        boolean = True
    elif s in {"f", "false", "0", "n", "no"}:
        boolean = False
    else:
        exit("Error: invalid argument specified for boolean flag: %s" % s)

    return boolean


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", stringAsBool)  # add type keyword to registries
    parser.add_argument(
        "--inputSequences",
        type=str,
        required=True,
        help="File path of FASTQ or FASTA sequence file containing sequences for assembly"
    )
    parser.add_argument(
        "--outputDir",
        type=str,
        default="./output/",
        required=False,
        help="Desired output directory path (will be created during run time if doesn't exist)"
    )
    parser.add_argument(
        "--minReadLength",
        type=int,
        # default=1000,
        required=False,
        help="The minimum read length. Reads shorter than this are skipped on input."
    )
    parser.add_argument(
        "--k",
        type=int,
        # default=10,
        required=False,
        help="The length of the k-mers used as markers. \n \
              When useRunLengthReads is False, suggested choice is k = 8. \n \
              When useRunLengthReads is True, suggested choice is k = 10."
    )
    parser.add_argument(
        "--probability",
        type=float,
        # default=0.1,
        required=False,
        help="The probability that a k-mer is a marker. \n \
              This is approximately equal to the fraction\n \
              of k-mers that will be used as markers."
    )
    parser.add_argument(
        "--m",
        type=int,
        # default=4,
        required=False,
        help="The number of consecutive markers that define a MinHash feature."
    )
    parser.add_argument(
        "--minHashIterationCount",
        type=int,
        # default=100,
        required=False,
        help="The number of MinHash iterations."
    )
    parser.add_argument(
        "--maxBucketSize",
        type=int,
        # default=30,
        required=False,
        help="The maximum bucket size to be used by the MinHash algoritm. \n \
              Buckets larger than this are ignored."
    )
    parser.add_argument(
        "--minFrequency",
        type=int,
        # default=1,
        required=False,
        help="The minimum number of times a pair of oriented reads \n \
              is found by the MinHash algorithm for the pair to \n \
              generate an overlap."
    )
    parser.add_argument(
        "--maxSkip",
        type=int,
        # default=30,
        required=False,
        help="The maximum number of markers that an alignment is allowed\n \
              to skip on either of the oriented reads being aligned."
    )
    parser.add_argument(
        "--maxMarkerFrequency",
        type=int,
        # default=10,
        required=False,
        help="Marker frequency threshold. \n \
              When computing an alignment between two oriented reads, \n \
              marker kmers that appear more than this number of times \n \
              in either of the two oriented reads are discarded \n \
              (in both oriented reads)."
    )
    parser.add_argument(
        "--minAlignedMarkerCount",
        type=int,
        # default=100,
        required=False,
        help="The minimum number of aligned markers in an alignment \n \
              in order for the alignment to be considered good and usable."
    )
    parser.add_argument(
        "--maxTrim",
        type=int,
        # default=30,
        required=False,
        help="The maximum number of trim markers tolerated at the \n \
              beginning and end of an alignment. There can be \n \
              up this number of markers between the first/last aligned marker \n \
              and the beginning/end of either oriented read \n \
              for an alignment to be considered good and usable."
    )
    parser.add_argument(
        "--minComponentSize",
        type=int,
        # default=100,
        required=False,
        help="The minimum size (number of oriented reads) of \n \
              a connected component to be kept."
    )
    parser.add_argument(
        "--maxChimericReadDistance",
        type=int,
        # default=2,
        required=False,
        help="Argument maxChimericReadDistance for flagChimericReads."
    )
    parser.add_argument(
        "--minCoverage",
        type=int,
        # default=10,
        required=False,
        help="The minimum and maximum coverage (number of markers) \n \
              for a vertex of the marker graph. \n \
              Vertices with coverage outside this range are collapsed \n \
              away and not generated by computeMarkerGraphVertices."
    )
    parser.add_argument(
        "--maxCoverage",
        type=int,
        # default=100,
        required=False,
        help="The minimum and maximum coverage (number of markers) \n \
              for a vertex of the marker graph. \n \
              Vertices with coverage outside this range are collapsed \n \
              away and not generated by computeMarkerGraphVertices."
    )
    parser.add_argument(
        "--lowCoverageThreshold",
        type=int,
        # default=1,
        required=False,
        help="Parameters for flagMarkerGraphWeakEdges."
    )
    parser.add_argument(
        "--highCoverageThreshold",
        type=int,
        # default=1000,
        required=False,
        help="Parameters for flagMarkerGraphWeakEdges."
    )
    parser.add_argument(
        "--maxDistance",
        type=int,
        # default=30,
        required=False,
        help="Parameters for flagMarkerGraphWeakEdges."
    )
    parser.add_argument(
        "--pruneIterationCount",
        type=int,
        # default=6,
        required=False,
        help="Number of iterations for pruneMarkerGraphStrongSubgraph."
    )
    parser.add_argument(
        "--shortCycleLengthThreshold",
        type=int,
        # default=30,
        required=False,
        help="Maximum length (number of edges) for a short cycle to be removed."
    )
    parser.add_argument(
        "--bubbleLengthThreshold",
        type=int,
        # default=1000,
        required=False,
        help="Maximum length (number of edges) for a bubble to be removed."
    )
    parser.add_argument(
        "--superBubbleLengthThreshold",
        type=int,
        # default=100,
        required=False,
        help="Length threshold (number of edges) used for superbubble removal."
    )
    parser.add_argument(
        "--markerGraphEdgeLengthThresholdForConsensus",
        type=int,
        # default=100,
        required=False,
        help="Used during sequence assembly."
    )
    parser.add_argument(
        "--consensusCaller",
        type=str,
        required=False,
        choices=["Simple","SimpleBayesian"],
        help="Whether to use Bayesian inference on read lengths during consensus calling"
    )
    parser.add_argument(
        "--useMarginPhase",
        type="bool",
        # default=True,
        required=False,
        help="Use margin polisher during consensus. \n\n \
              Any case insensitive variant of the following is accepted: \n \
              t, true, 1, y, yes, f, false, 0, n, no"
    )

    args = parser.parse_args()

    main(readsSequencePath=args.inputSequences, 
         outputParentDirectory=args.outputDir,
         args=args)
