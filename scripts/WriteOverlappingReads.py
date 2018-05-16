#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = """
This will write in fasta format all the oriented
reads tghat overlap a given oriented read.
Invoke with two arguments: read id and strand.
"""
if not len(sys.argv)==3:
    print(helpMessage)
    exit(1)
    
readId = int(sys.argv[1])
strand= int(sys.argv[2])

   

a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessOverlaps()
a.writeOverlappingReads(readId=readId, strand=strand)

