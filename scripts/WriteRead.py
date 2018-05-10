#!/usr/bin/python3

import Nanopore2

readId = int(input('Enter a ReadId: '))
fileName = str(readId) + '.fasta'

a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.writeRead(readId=readId, fileName=fileName)

