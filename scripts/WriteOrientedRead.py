#!/usr/bin/python3

import Nanopore2
import sys

helpMessage = """
Invoke with two arguments, the id and strand of the read to be written out.
"""

if not len(sys.argv) == 3:
    print(helpMessage)
    exit(1)

readId = int(sys.argv[1])
strand = int(sys.argv[2])
if not (strand==0 or strand==1):
    print('Invalid strand')
    print(helpMessage)
    exit(1)


fileName = str(readId) + '-' + str(strand) + '.fasta'

a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.writeOrientedRead(readId=readId, strand=strand, fileName=fileName)

