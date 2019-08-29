#!/usr/bin/python3

# Import what we need.
import argparse
import shasta

# Get the arguments.
parser = argparse.ArgumentParser(description="Compute the phasing similarity between two oriented reads.")
parser.add_argument("readId0", type=int, help="The read id of the first oriented read.")
parser.add_argument("strand0", type=int, help="The strand (0 or 1) of the first oriented read.")
parser.add_argument("readId1", type=int, help="The read id of the second oriented read.")
parser.add_argument("strand1", type=int, help="The strand (0 or 1) of the second oriented read.")
arguments = parser.parse_args()


# Create the Assembler.
a = shasta.Assembler()

# Access what we need.
a.accessPhasingGraph();

# Do the work.
print(a.computePhasingSimilarity(
    arguments.readId0, arguments.strand0,
    arguments.readId1, arguments.strand1))





