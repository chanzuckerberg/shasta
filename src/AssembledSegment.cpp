#include "AssembledSegment.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



void AssembledSegment::clear()
{
    assemblyGraphEdgeId = AssemblyGraph::invalidEdgeId;
    vertexCount = 0;
    edgeCount = 0;

    runLengthSequence.clear();
    repeatCounts.clear();
    vertexIds.clear();
    vertexCoverage.clear();
}
