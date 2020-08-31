#!/usr/bin/python3

# This runs a standard sequence of assembly steps, assuming
# that the read graph was already created.
# Only useful for debugging.

import runpy

runpy.run_module("CreateMarkerGraphVertices")
runpy.run_module("FindMarkerGraphReverseComplementVertices")
runpy.run_module("CreateMarkerGraphEdges")
runpy.run_module("FindMarkerGraphReverseComplementEdges")
runpy.run_module("TransitiveReduction")
runpy.run_module("PruneMarkerGraphStrongSubgraph")
runpy.run_module("SimplifyMarkerGraph")
runpy.run_module("CreateAssemblyGraphEdges")
runpy.run_module("CreateAssemblyGraphVertices")
runpy.run_module("AssembleMarkerGraphVertices")
runpy.run_module("AssembleMarkerGraphEdges")
runpy.run_module("Assemble")
runpy.run_module("ComputeAssemblyStatistics")
runpy.run_module("WriteGfaBothStrands")

