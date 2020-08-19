#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessAlignmentData()
a.accessReadGraph()
a.removeReadGraphBridges()

