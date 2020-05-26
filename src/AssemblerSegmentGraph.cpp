#include "Assembler.hpp"
#include "SegmentGraph.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

#include <map>

// The segment graph is used to detangle and analyze reachability
// in the assembly graph.

// In the segment graph, each vertex corresponds to an assembly
// graph edge (segment). But not every assembly graph edge corresponds to
// a segment graph vertex.
// An edge is created between two vertices if there are one
// or more oriented reads that encounter the two corresponding key
// assembly graph edges consecutively and in that order.


void Assembler::createSegmentGraph()
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
        " vertices and " << assemblyGraph.edges.size() << " edges." << endl;




    // Each oriented read is guaranteed to correspond to a path in the marker graph.
    // But because of transitive reduction, that path does not necessarily correspond
    // to a path in the assembly graph.
    // Here we find, for each oriented read, the sequence of assembly graph edges
    // encountered along its marker graph path.
    // This vector is indexed by OrientedReadId::getValue().
    vector< vector<AssemblyGraph::EdgeId> > orientedReadAssemblyGraphEdges(2*readCount());
    vector<MarkerGraph::EdgeId> markerGraphPath;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            vector<AssemblyGraph::EdgeId>& assemblyGraphEdges =
                orientedReadAssemblyGraphEdges[orientedReadId.getValue()];

            // Find the marker graph path of this oriented read.
            computeOrientedReadMarkerGraphPath(
                orientedReadId,
                0, uint32_t(markers.size(orientedReadId.getValue())-1),
                markerGraphPath);

            // Loop over the path.
            AssemblyGraph::EdgeId previousAssemblyGraphEdgeId =
                std::numeric_limits<AssemblyGraph::EdgeId>::max();
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

                // There is only one assembly graph edge.
                const AssemblyGraph::EdgeId assemblyGraphEdgeId = v.front().first;

                // If same as the previous, slip.
                if(assemblyGraphEdgeId == previousAssemblyGraphEdgeId) {
                    continue;
                }

                // This is the next key assembly graph edge encountered
                // by this oriented read along its marker graph path. Store it.
                assemblyGraphEdges.push_back(assemblyGraphEdgeId);
                previousAssemblyGraphEdgeId = assemblyGraphEdgeId;
            }
        }
    }



    // Write a csv file with the sequence of assembly graph edges
    // for each oriented read.
    ofstream csv("OrientedReadSegments.csv");
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            csv << orientedReadId << ",";

            vector<AssemblyGraph::EdgeId>& keyAssemblyGraphEdges =
                orientedReadAssemblyGraphEdges[orientedReadId.getValue()];
            for(const AssemblyGraph::EdgeId edgeId: keyAssemblyGraphEdges) {
                csv << edgeId << ",";
            }
            csv << "\n";
        }
    }



    // Initially, we will generate segment graph vertices only for
    // assembly graph edges (segments) v0->v1 such that in-degree(v0)<2 and
    // out_degree(v)<2, that is, there is no uncertainty on what precedes
    // and follows the segment. In some other places in the code, these
    // are called "key" edges of the assembly graph or "key" segments.
    vector<bool> isSegmentIncluded(assemblyGraph.edges.size());
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        const AssemblyGraph::VertexId v0 = edge.source;
        const AssemblyGraph::VertexId v1 = edge.target;
        const uint64_t inDegree0 = assemblyGraph.edgesByTarget.size(v0);
        const uint64_t outDegree1 = assemblyGraph.edgesBySource.size(v1);
        isSegmentIncluded[edgeId] = (inDegree0 < 2) and (outDegree1 < 2);
    }



    // Main iteration. At each iteration, some vertices disappear.
    for(uint64_t iteration=0; ; iteration++) {
        cout << "Iteration " << iteration << " begins." << endl;

        // Create SegmentGraph vertices.
        // We generate a vertex for each assembly graph edge (segment)
        // for which isSegmentIncluded is true.
        segmentGraphPointer = make_shared<SegmentGraph>();
        SegmentGraph& segmentGraph = *segmentGraphPointer;
        std::map<AssemblyGraph::EdgeId, SegmentGraph::vertex_descriptor> vertexMap;
        for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
            if(isSegmentIncluded[edgeId]) {
                SegmentGraph::vertex_descriptor v = add_vertex(SegmentGraphVertex(edgeId), segmentGraph);
                vertexMap.insert(make_pair(edgeId, v));
            }
        }

        // Create segment graph edges.
        // Follow the sequence of assembly graph edges
        // for each oriented read.
        vector<AssemblyGraph::EdgeId> v;
        for(ReadId readId=0; readId<readCount(); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);

                // The sequence of assembly graph edges encountered by this oriented read.
                vector<AssemblyGraph::EdgeId>& assemblyGraphEdges =
                    orientedReadAssemblyGraphEdges[orientedReadId.getValue()];

                // Gather into working vector v the assembly graph edges
                // of this read for which isSegmentIncluded is true.
                v.clear();
                for(const AssemblyGraph::EdgeId edgeId: assemblyGraphEdges) {
                    if(isSegmentIncluded[edgeId]) {
                        v.push_back(edgeId);
                    }
                }

                // Now we can generate the edges.
                for(uint64_t i=1; i<v.size(); i++) {
                    const AssemblyGraph::EdgeId edgeId0 = v[i-1];
                    const AssemblyGraph::EdgeId edgeId1 = v[i];

                    // Find the corresponding vertices of the segment graph.
                    const auto it0 = vertexMap.find(edgeId0);
                    SHASTA_ASSERT(it0 != vertexMap.end());
                    const SegmentGraph::vertex_descriptor v0 = it0->second;
                    const auto it1 = vertexMap.find(edgeId1);
                    SHASTA_ASSERT(it1 != vertexMap.end());
                    const SegmentGraph::vertex_descriptor v1 = it1->second;

                    // Add this segment graph edge, if we don't already have it.
                    SegmentGraph::edge_descriptor e;
                    bool edgeExists;
                    tie(e, edgeExists) = edge(v0, v1, segmentGraph);
                    if(not edgeExists) {
                        tie(e, ignore) = add_edge(v0, v1, segmentGraph);
                    }
                    segmentGraph[e].coverage++;
                }
            }
        }

        // Remove low coverage edges.
        const uint64_t minCoverage = 2; // EXPOSE WHEN CODE STABILIZES. *********
        segmentGraph.removeLowCoverageEdges(minCoverage);

        // Approximate transitive reduction.
        const uint64_t transitiveReductionMaxDistance = 4;   // EXPOSE WHEN CODE STABILIZES. *********
        segmentGraph.transitiveReduction(transitiveReductionMaxDistance);

        segmentGraph.writeGraphviz("SegmentGraph-" + to_string(iteration) + ".dot");
        cout << "The segment graph has " << num_vertices(segmentGraph) <<
            " vertices and " << num_edges(segmentGraph) << " edges." << endl;


        // Keep only segments for which both in-degree and out-degree is less than 2.
        uint64_t removedCount = 0;
        BGL_FORALL_VERTICES(v, segmentGraph, SegmentGraph) {
            if((in_degree(v, segmentGraph)>1) or (out_degree(v, segmentGraph)>1)) {
                isSegmentIncluded[segmentGraph[v].assemblyGraphEdgeId] = false;
                removedCount++;
            }
        }
        cout << "Removed " << removedCount << " segments." << endl;
        if(removedCount > 0) {
            continue;
        }

        // When we get here, all vertices of the segment graph have in-degree and
        // out-degree less than 2. This means that the segment graph does not
        // contain any branches. Therefore, each connected component
        // of the segment graph is a linear chain, possibly circular, and
        // possibly consisting of a single vertex.
        segmentGraph.findChains();
        return;
    }
}



void Assembler::colorGfaBySegmentGraphChain(uint64_t chainId) const
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    const SegmentGraph& segmentGraph = *segmentGraphPointer;

    // Find the chain each segment belongs to, if any.
    const uint64_t noChain = std::numeric_limits<uint64_t>::max();
    vector<uint64_t> chainTable(assemblyGraph.edges.size(), noChain);
    for(uint64_t chainId=0; chainId<segmentGraph.chains.size(); chainId++) {
        const SegmentGraph::Chain& chain = segmentGraph.chains[chainId];
        for(const SegmentGraph::vertex_descriptor v: chain.vertices) {
            const AssemblyGraph::EdgeId edgeId = segmentGraph[v].assemblyGraphEdgeId;
            chainTable[edgeId] = chainId;
        }
    }

    ofstream csv("AssemblyGraph-BothStrands-Color.csv");
    csv << "AssemblyGraphEdgeId,ChainId,Color\n";
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<chainTable.size(); edgeId++) {
        const auto c = chainTable[edgeId];
        csv << edgeId << ",";
        if(c != noChain) {
            csv << c;
        }
        csv << ",";
        if(c == chainId) {
            csv << "Blue";
        } else {
            csv << "Grey";
        }
        csv << "\n";
    }

}
