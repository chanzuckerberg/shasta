#!/usr/bin/python3

import shasta
import GetConfig
import sys

# Read the config file.
config = GetConfig.getConfig()

helpMessage = "Invoke with one argument, the name of the Fasta file."

if not len(sys.argv)==2:
    print(helpMessage)
    exit(1)
    
fileName = sys.argv[1]

a = shasta.Assembler()
a.addReadsFromFasta(
    fileName = fileName, 
    minReadLength = int(config['Reads']['minReadLength']))

