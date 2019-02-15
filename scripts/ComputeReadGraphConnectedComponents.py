#!/usr/bin/python3

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessReadsReadOnly()
a.accessAlignmentData()
a.accessReadGraph()
a.accessReadFlags(readWriteAccess = True)
a.computeReadGraphConnectedComponents(
    int(config['ReadGraph']['minComponentSize']))

