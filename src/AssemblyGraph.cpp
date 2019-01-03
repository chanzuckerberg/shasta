#include "AssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Close and remove all open data.
void AssemblyGraph::remove()
{
    if(vertices.isOpen()) {
        vertices.remove();
    }

    if(edges.isOpen) {
        edges.remove();
    }
    if(edgesBySource.isOpen()) {
        edgesBySource.remove();
    }
    if(edgesByTarget.isOpen()) {
        edgesByTarget.remove();
    }
    if(markerToAssemblyTable.isOpen) {
        markerToAssemblyTable.remove();
    }
    if(sequences.isOpen()) {
        sequences.remove();
    }
    if(repeatCounts.isOpen()) {
        repeatCounts.remove();
    }
}
