#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessReadsReadOnly()
a.accessKmers()
a.findMarkers()

