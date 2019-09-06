#!/usr/bin/python3

# Import what we need.
import argparse
import shasta

# Get the arguments.
parser = argparse.ArgumentParser(description="Compute the phasing similarity between two assembly graph edges.")
parser.add_argument("edgeId0", type=int, help="The id of the first assembly graph edge.")
parser.add_argument("edgeId1", type=int, help="The id of the second assembly graph edge.")
arguments = parser.parse_args()


# Create the Assembler.
a = shasta.Assembler()

# Access what we need.
a.accessPhasingGraph();

# Do the work.
print('Jaccard similarity %f' % a.computePhasingSimilarity(arguments.edgeId0, arguments.edgeId1))
print('Number of common reads %i' % a.countCommonInternalOrientedReads(arguments.edgeId0, arguments.edgeId1))






