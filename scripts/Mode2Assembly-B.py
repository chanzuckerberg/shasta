#!/usr/bin/python3

"""

This run the final portion of Mode 2 assembly.

"""

import ast
import shasta
import GetConfig

config = GetConfig.getConfig()

shasta.openPerformanceLog('Mode2Assembly-B.log')

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()

mode2Options = shasta.Mode2AssemblyOptions();
mode2Options.suppressGfaOutput      = ast.literal_eval(config['Assembly']['mode2.suppressGfaOutput'])
mode2Options.suppressFastaOutput    = ast.literal_eval(config['Assembly']['mode2.suppressFastaOutput'])
mode2Options.suppressDetailedOutput = ast.literal_eval(config['Assembly']['mode2.suppressDetailedOutput'])
mode2Options.suppressPhasedOutput   = ast.literal_eval(config['Assembly']['mode2.suppressPhasedOutput'])
mode2Options.suppressHaploidOutput  = ast.literal_eval(config['Assembly']['mode2.suppressHaploidOutput'])

a.createAssemblyGraph2(
    superbubbleRemovalEdgeLengthThreshold = int(config['Assembly']['superbubbleRemoval.edgeLengthThreshold']),
    pruneLength = int(config['Assembly']['pruneLength']),
    phasingMinReadCount = int(config['Assembly']['phasing.minReadCount']),
    mode2Options = mode2Options 
    )




