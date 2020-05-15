#include "Assembler.hpp"
using namespace shasta;

// The segment graph is used to detangle and determine reachability
// in the assembly graph.



void Assembler::createSegmentGraph()
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;

    // Figure out which segments are essential.
    // A segment (assembly graph edge) v0->v1 is a key segment if in-degree(v0)<2 and
    // out_degree(v)<2, that is, there is no uncertainty on what preceeds
    // and follows the segment.
    vector<bool> isKeyEdge(assemblyGraph.edges.size());
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        const AssemblyGraph::VertexId v0 = edge.source;
        const AssemblyGraph::VertexId v1 = edge.target;
        const uint64_t inDegree0 = assemblyGraph.edgesByTarget.size(v0);
        const uint64_t outDegree1 = assemblyGraph.edgesBySource.size(v1);
        isKeyEdge[edgeId] = (inDegree0 < 2) and (outDegree1 < 2);
    }



    // Loop over all oriented reads.
    vector<MarkerGraph::EdgeId> markerGraphPath;
    std::map<pair<MarkerGraph::EdgeId, MarkerGraph::EdgeId>, uint64_t> segmentGraphEdges;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);

            // Find the marker graph path of this oriented read.
            computeOrientedReadMarkerGraphPath(
                orientedReadId,
                0, uint32_t(markers.size(orientedReadId.getValue())-1),
                markerGraphPath);

            // Loop over the path, looking for key assembly graph segments (edges).
            AssemblyGraph::EdgeId previousAssemblyGraphEdgeId = std::numeric_limits<AssemblyGraph::EdgeId>::max();
            for(const MarkerGraph::EdgeId markerGraphEdgeId: markerGraphPath) {

                // Get the corresponding assembly graph edges.
                const span<const pair<AssemblyGraph::EdgeId, uint32_t> > v =
                    assemblyGraph.markerToAssemblyTable[markerGraphEdgeId];

                // If no assembly graph edges, skip.
                if(v.size() == 0) {
                    continue;
                }

                // If detangling was used, there can be more than one,
                // and we don't want this here.
                SHASTA_ASSERT(v.size() == 1);

                // Of, there is only one.
                const AssemblyGraph::EdgeId assemblyGraphEdgeId = v.front().first;

                // If not a key edge, skip.
                if(not isKeyEdge[assemblyGraphEdgeId]) {
                    continue;
                }

                // If same as the previous, slip.
                if(assemblyGraphEdgeId == previousAssemblyGraphEdgeId) {
                    continue;
                }

                if(previousAssemblyGraphEdgeId != std::numeric_limits<AssemblyGraph::EdgeId>::max()) {
                    // We found a transition between key segments.
                    // This generates an edge in the segment graph.
                    const auto key = make_pair(previousAssemblyGraphEdgeId, assemblyGraphEdgeId);
                    const auto it = segmentGraphEdges.find(key);
                    if(it == segmentGraphEdges.end()) {
                        segmentGraphEdges.insert(
                            make_pair(make_pair(previousAssemblyGraphEdgeId, assemblyGraphEdgeId), 1));
                    } else {
                        ++it->second;
                    }
                }

                previousAssemblyGraphEdgeId = assemblyGraphEdgeId;
            }

        }
    }

    ofstream graphOut("SegmentGraph.dot");
    graphOut << "digraph G {\n";
    for(const auto& p: segmentGraphEdges) {
        if(p.second > 0) {
            graphOut << p.first.first << "->" << p.first.second <<
                " [penwidth=\"" << 0.2*double(p.second) << "\"];\n";
        }
    }
    graphOut << "}\n";
}
