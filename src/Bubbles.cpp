#include "Bubbles.hpp"
using namespace shasta;


Bubbles::Bubbles(
    const AssemblyGraph& assemblyGraph) :
    assemblyGraph(assemblyGraph)
{
    findBubbles();
    cout << "Found " << bubbles.size() << " bubbles." << endl;
}



// For now we only consider diploid bubbles, defined using the
// following strict criteria:
// - Source vertex v0 has out-degree 2.
// - Target vertex v1 has in-degree 2.
// - There are two parallel edges eA and eB, both v0->v1.

void Bubbles::findBubbles()
{
    // Loop over possible choices for v0 (source vertex of the bubble).
    const AssemblyGraph::VertexId vertexCount = assemblyGraph.vertices.size();
    for(AssemblyGraph::VertexId v0=0; v0<vertexCount; v0++) {

        // Check the out-degree.
        const auto outEdges0 = assemblyGraph.edgesBySource[v0];
        if(outEdges0.size() != 2) {
            continue;
        }

        const AssemblyGraph::EdgeId eA = outEdges0[0];
        const AssemblyGraph::EdgeId eB = outEdges0[1];
        const AssemblyGraph::Edge& edgeA = assemblyGraph.edges[eA];
        const AssemblyGraph::Edge& edgeB = assemblyGraph.edges[eB];

        // Check the target vertex.
        const AssemblyGraph::VertexId v1 = edgeA.target;
        if(edgeB.target != v1) {
            continue;
        }

        // We have a diploid bubble.
        bubbles.push_back(Bubble());

    }
}
