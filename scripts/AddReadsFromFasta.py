#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = "Invoke with the one argument, the name of the Fasta file."

if not len(sys.argv)==2:
    print(helpmessage)
    exit(1)
    
fileName = sys.argv[1]

a = Nanopore2.Assembler()
a.accessReadsReadWrite();
a.accessReadNamesReadWrite();
a.addReadsFromFasta(fileName=fileName)

