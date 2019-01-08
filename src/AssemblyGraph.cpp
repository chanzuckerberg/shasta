#include "AssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

#include "fstream.hpp"

// Close and remove all open data.
void AssemblyGraph::remove()
{
    if(vertices.isOpen) {
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

    if(edgeLists.isOpen()) {
        edgeLists.remove();
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


// Basic Graphviz output of the global assembly graph.
void AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream graphOut(fileName);
    graphOut << "digraph AssemblyGraph {\n";
    for(const Edge& edge: edges) {
        graphOut << edge.source << "->" << edge.target << ";\n";
    }
    graphOut << "}\n";
}
