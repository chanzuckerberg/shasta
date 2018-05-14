#!/usr/bin/python3

import Nanopore2
import sys

   

a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.accessOverlaps()
a.computeOverlapGraphComponents(
    minFrequency = 5,
    minComponentSize = 100
    )

