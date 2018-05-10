#!/usr/bin/python3

import Nanopore2

a = Nanopore2.Assembler()
a.accessReadsReadOnly()
a.accessKmers()
a.findMarkers()

