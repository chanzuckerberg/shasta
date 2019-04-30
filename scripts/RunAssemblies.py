#!/usr/bin/python3

import glob
import os

"""

This script can be used to run multiple assemblies.
It assumes that the Shasta executable is in the
same directory that contains this script
(typically, that will be the shasta-install/bin directory).

The current directory must contain a file called
samples that lists the sample names to be assembled, one per line.
Blank lines and lines beginning wit "#" are skipped.

The current directory must also contain all of
the input fasta files for all the samples.
The input files for a sample named xyz must be named

xyz*.fasta

"""

# Locate the Shasta executable. It is assumed to be in the
# same directory that contains this script.
shastaExecutable = os.path.dirname(os.path.realpath(__file__)) + '/shasta'
if not os.path.exists(shastaExecutable):
    raise Exception('Shasta executable not found at %s' % shastaExecutable)
    

# Read in the sample names.
samples = []
for line in open('samples', 'r'):

    # Remove the line end and leading/trailing white space.
    line = line.rstrip()
    line = line.lstrip()
    
    # Remove leading and trailing space.
    
    # Skip empty lines and lines beginning with "#".
    if not line:
        continue
    if line[0] == '#':
        continue

    # This is a sample. Store it.
    sample = line
    samples.append(sample) 
    
print('Found %i samples:' % len(samples))
for sample in samples:
    print(sample)  
    


# Main loop over samples.
for sample in samples:
    
    # Get the input file names for this sample.
    inputFileNames = glob.glob(sample + '*.fasta')
    print('Working on %s with %i input files:' % (sample, len(inputFileNames))) 
    for inputFileName in inputFileNames:
        print(inputFileName)
        
    # Run this sample.
            

