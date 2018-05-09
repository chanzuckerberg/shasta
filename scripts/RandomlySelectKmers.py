#!/usr/bin/python3

import Nanopore2

a = Nanopore2.Assembler()
a.randomlySelectKmers(k=8, probability=0.1)
a.writeKmers()
