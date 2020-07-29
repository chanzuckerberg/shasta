#!/usr/bin/python3

import shasta
import sys




a = shasta.Assembler()
a.getReads().writeReads(fileName='Reads.fasta')

