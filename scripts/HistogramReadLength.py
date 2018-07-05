#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessReadsReadOnly();
a.accessReadNamesReadOnly();
a.histogramReadLength(fileName="ReadLengthHistogram.csv")

