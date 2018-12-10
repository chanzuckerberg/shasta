#include "AssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Write the assembly graph in GFA 1.0 format defined here:
// https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
void AssemblyGraph::writeGfa1(ostream& gfa)
{
    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment for each vertex.
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
        const auto sequence = sequences[vertexId];
        gfa << "S\t" << vertexId << "\t";
        gfa << sequence << "\n";
    }

    // Write a link for each edge.
    for(const Edge& edge: edges) {
        gfa << "L\t" <<
            edge.source << "\t" <<
            "+\t" <<
            edge.target << "\t" <<
            "+\t" <<
            "*\n";
    }
}
