#!/usr/bin/python3

import configparser
import sys
import getopt
import json
import csv

helpMessage = """
Usage:
    GenerateFeedback.py --assemblyDirectory /path/to/assemblyDirectory

Note: The goal of this script is to provide feedback on a de novo assembly done using the
Shasta long read assembler.

Each genome is different and reads differ in quality. So it is likely that you will need
to repeat the `assembly -> feedback -> assembly` process a few times.

If you're not happy with the assembly after a few (3-4) tries, please file a Github issue
and we can help.
"""

def usage():
    print(helpMessage)
    return


def loadAssemblySummary(assemblyDirPath):
    assemblySummary = {}
    assemblySummaryFile = '{}/AssemblySummary.json'.format(assemblyDirPath)
    with open(assemblySummaryFile) as jsonFile:
        assemblySummary = json.load(jsonFile)

    return assemblySummary


def analyze(assemblyDirPath, genomeSize):
    assemblySummary = loadAssemblySummary(assemblyDirPath)
    readsUsedInAssembly = assemblySummary['Reads used in this assembly']
    numberOfReads = readsUsedInAssembly['Number of reads']

    readGraph = assemblySummary['Read graph']
    isolatedReadsFraction = float(readGraph['Isolated reads fraction']['Reads'])

    alignments = assemblySummary['Alignments']
    numberOfAlignmentCandidates = alignments['Number of alignment candidates found by the LowHash algorithm']
    numberOfGoodAlignments = alignments['Number of good alignments']

    assembledSegments = assemblySummary['Assembled segments']
    totalAssembledLength = assembledSegments['Total assembled segment length']
    longestSegmentLength = assembledSegments['Longest assembled segment length']
    segmentsN50 = assembledSegments['Assembled segments N50']

    print()
    print('Number of reads used = {:,}'.format(numberOfReads))
    print('Isolated reads fraction = {:,.2f}'.format(isolatedReadsFraction))
    print('Number of alignment candidates = {:,}'.format(numberOfAlignmentCandidates))
    print('Number of good alignments = {:,}'.format(numberOfGoodAlignments))
    print()
    print('Genome fraction assembled = {:,.2f} %'.format(totalAssembledLength * 100 / (genomeSize * 1024 * 1024)))
    print('Longest assembled segment length = {:,}'.format(longestSegmentLength))
    print('Assembled segments N50 = {:,}'.format(segmentsN50))
    print()

    avgCandidatesPerRead = numberOfAlignmentCandidates / numberOfReads

    config = getConfig(assemblyDirPath)
    minHashConfig = config['MinHash']
    
    print('Feedback, if any:')

    if (avgCandidatesPerRead < 20):
        print('MinHash phase did not generate enough alignment candidates.')
        print('Try the following in order:')
        print('  (Suggestion) Increase `MinHash.minHashIterationCount` by 10, up to a maximum of 100.')
        if (int(minHashConfig['m']) == 4):
            print('  (Suggestion) Decrease `MinHash.m` to 3.')
    
    else:
        # Enough promising candidate pairs were generated
        avgGoodAlignmentsPerRead = numberOfGoodAlignments / numberOfReads
        if (avgGoodAlignmentsPerRead < 5 or isolatedReadsFraction > 0.5):
            # ... but not enough candidates met the bar of what is considered a good alignment.
            msg = (
                'Not enough good alignments were generated per read. '
                'Try relaxing the definition of what makes a good alignment.'
            )
            print(msg)
            print('Try the following in order:')
            print('  (Suggestion) Decrease `Align.minAlignedFraction` by 0.05, up to a minimum of 0.2.')
            print('  (Suggestion) Decrease `Align.minAlignedMarkerCount` by 20, up to a minimum of 200.')
            print('  (Suggestion) Increase `Align.maxSkip` & `Align.maxDrift` by 10, to allow for larger gaps in alignments.')
        

def getConfig(assemblyDirPath):
    config = configparser.ConfigParser()
    configFilePath = '{}/shasta.conf'.format(assemblyDirPath)
    if not config.read(configFilePath):
        raise Exception('Error reading config file {}.'.format(configFilePath))
    return config


def main(argv):
    assemblyDirPath = ""
    try:
      opts, args = getopt.getopt(argv,"", ["assemblyDirectory="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
      if opt in ("--assemblyDirectory"):
         assemblyDirPath = arg

    if assemblyDirPath == "":
        usage()
        exit(2)

    print()
    genomeSize = int(input('Approximate genome size in megabasis (Mbp): '))
    analyze(assemblyDirPath, genomeSize)
    return

if __name__ == '__main__':
    main(sys.argv[1:])
