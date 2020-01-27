#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessDirectedReadGraphReadWrite()
a.accessConflictReadGraph()
a.markDirectedReadGraphConflictEdges1()


