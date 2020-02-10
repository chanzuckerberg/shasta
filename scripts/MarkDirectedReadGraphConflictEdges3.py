#!/usr/bin/python3

import shasta

a = shasta.Assembler()
a.accessDirectedReadGraphReadWrite()
a.accessConflictReadGraph()

# For testing use a fixed radius.
a.markDirectedReadGraphConflictEdges3(radius = 6)


