#!/usr/bin/python3

import shasta
import sys

helpMessage = """
Invoke with one argument, the id of the read to be written out.
"""

if not len(sys.argv) == 2:
    print(helpMessage)
    exit(1)

readId = int(sys.argv[1])
fileName = str(readId) + '.fasta'

a = shasta.Assembler()
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.writeRead(readId=readId, fileName=fileName)

