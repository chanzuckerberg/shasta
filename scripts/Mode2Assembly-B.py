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



# Fill in the Mode2Assemblyoptions.
mode2Options = shasta.Mode2AssemblyOptions();

mode2Options.strongBranchThreshold = int(config['Assembly']['mode2.strongBranchThreshold'])
mode2Options.epsilon = float(config['Assembly']['mode2.epsilon'])

mode2Options.minConcordantReadCountForBubbleRemoval = int(config['Assembly']['mode2.bubbleRemoval.minConcordantReadCount'])
mode2Options.maxDiscordantReadCountForBubbleRemoval = int(config['Assembly']['mode2.bubbleRemoval.maxDiscordantReadCount'])
mode2Options.minLogPForBubbleRemoval = float(config['Assembly']['mode2.bubbleRemoval.minlogP'])

mode2Options.componentSizeThresholdForBubbleRemoval = int(config['Assembly']['mode2.bubbleRemoval.componentSizeThreshold'])
mode2Options.minConcordantReadCountForPhasing = int(config['Assembly']['mode2.phasing.minConcordantReadCount'])
mode2Options.maxDiscordantReadCountForPhasing = int(config['Assembly']['mode2.phasing.maxDiscordantReadCount'])
mode2Options.minLogPForPhasing = float(config['Assembly']['mode2.phasing.minlogP'])

mode2Options.maxSuperbubbleSize = int(config['Assembly']['mode2.superbubble.maxSize'])
mode2Options.maxSuperbubbleChunkSize = int(config['Assembly']['mode2.superbubble.maxChunkSize'])
mode2Options.maxSuperbubbleChunkPathCount = int(config['Assembly']['mode2.superbubble.maxChunkPathCount'])
mode2Options.superbubbleEdgeLengthThreshold = int(config['Assembly']['mode2.superbubble.edgeLengthThreshold'])

mode2Options.suppressGfaOutput      = ast.literal_eval(config['Assembly']['mode2.suppressGfaOutput'])
mode2Options.suppressFastaOutput    = ast.literal_eval(config['Assembly']['mode2.suppressFastaOutput'])
mode2Options.suppressDetailedOutput = ast.literal_eval(config['Assembly']['mode2.suppressDetailedOutput'])
mode2Options.suppressPhasedOutput   = ast.literal_eval(config['Assembly']['mode2.suppressPhasedOutput'])
mode2Options.suppressHaploidOutput  = ast.literal_eval(config['Assembly']['mode2.suppressHaploidOutput'])




a.createAssemblyGraph2(
    pruneLength = int(config['Assembly']['pruneLength']),
    mode2Options = mode2Options 
    )




