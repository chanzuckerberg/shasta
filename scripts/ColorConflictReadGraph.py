#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessDirectedReadGraphReadOnly()
a.accessConflictReadGraph()
a.colorConflictReadGraph()


