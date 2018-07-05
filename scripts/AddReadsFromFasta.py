#!/usr/bin/python3

import Nanopore2
import Nanopore2GetConfig
import sys

# Read the config file.
config = Nanopore2GetConfig.getConfig()

helpMessage = "Invoke with the one argument, the name of the Fasta file."

if not len(sys.argv)==2:
    print(helpMessage)
    exit(1)
    
fileName = sys.argv[1]

a = Nanopore2.Assembler()
a.accessReadsReadWrite();
a.accessReadNamesReadWrite();
a.addReadsFromFasta(
    fileName = fileName, 
    minReadLength = int(config['Reads']['minReadLength']))

