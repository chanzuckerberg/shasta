#!/usr/bin/python3

import datetime
import glob
import os
import shutil

"""

This script can be used to run multiple assemblies.
It assumes that the Shasta executable is in the
same directory that contains this script
(typically, that will be the shasta-install/bin directory).

This script requires that the user can run sudo commands
without being prompted for a password.

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
    
    # Check that we continue to have sudo access.
    if not os.system('sudo --non-interactive true') == 0:
        raise Exception('Unable to run sudo.')
        
    # Get the input file names for this sample.
    inputFileNames = glob.glob(sample + '*.fasta')
    print(datetime.datetime.now())
    print('Working on %s with %i input files:' % (sample, len(inputFileNames))) 
    for inputFileName in inputFileNames:
        print(inputFileName, flush=True)
        
    # Clear Linux caches before running each sample.
    os.system('sudo sh -c "sync; echo 3 > /proc/sys/vm/drop_caches"')

    # Run this sample.
    command = ' '.join([
        shastaExecutable,
        '--memoryMode', 'filesystem',
        '--memoryBacking', '2M', 
        '--output', sample,
        '--input'] + inputFileNames) + ' 1>' + sample + '.stdout 2>&1'
    os.system(command)
    
    # Copy the log to the output directory.
    shutil.copy(sample + '.stdout', sample + '/' + sample + '.stdout')
    
    # Check that we continue to have sudo access.
    if not os.system('sudo --non-interactive true') == 0:
        raise Exception('Unable to run sudo.')

    # Cleanup the memory.
    command = ' '.join([
        shastaExecutable,
        '--command', 'cleanup',
        '--output', sample])
    os.system(command)
    
            

