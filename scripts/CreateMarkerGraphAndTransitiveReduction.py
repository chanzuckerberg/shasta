#!/usr/bin/python3

import runpy

runpy.run_module("CreateMarkerGraphVertices")
runpy.run_module("CreateMarkerGraphEdges")
runpy.run_module("FindMarkerGraphReverseComplementVertices")
runpy.run_module("FindMarkerGraphReverseComplementEdges")
runpy.run_module("TransitiveReduction")

