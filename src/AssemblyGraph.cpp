#include "AssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace ChanZuckerberg::shasta;

#include "fstream.hpp"

// Close and remove all open data.
void AssemblyGraph::remove()
{
    if(vertices.isOpen) {
        vertices.remove();
    }

    if(reverseComplementVertex.isOpen) {
    	reverseComplementVertex.remove();
    }

    if(edges.isOpen) {
        edges.remove();
    }

    if(reverseComplementEdge.isOpen) {
    	reverseComplementEdge.remove();
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

    // Write the vertices.
    // The label contains the corresponding marker graph vertex id.
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
        graphOut <<
            vertexId <<
            " [label=\"" <<
            vertexId << "\\n" << vertices[vertexId] <<
             "\"];\n";
    }

    // Write the edges.
    // The label contains the edge id and the number of maker graph edges
    // that correspond to this assembly graph edge.
    for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
        const Edge& edge = edges[edgeId];
        graphOut <<
            edge.source << "->" << edge.target <<
            " [label=\"" << edgeId << "\\n" <<
            edgeLists.size(edgeId) <<
            "\"];\n";
    }

    graphOut << "}\n";
}
