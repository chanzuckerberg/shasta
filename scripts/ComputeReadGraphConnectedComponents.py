#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessReadsReadOnly()
a.accessAlignmentData()
a.accessReadGraph()
a.accessChimericReadsFlags()
a.computeReadGraphConnectedComponents()

