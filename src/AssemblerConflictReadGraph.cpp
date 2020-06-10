#include "Assembler.hpp"
#include "approximateColoring.hpp"
#include "DynamicConflictReadGraph.hpp"
#include "makeBiconnected.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <map>
#include <set>
#include "tuple.hpp"



void Assembler::createConflictReadGraph(
    uint64_t threadCount,
    uint32_t maxOffsetSigma,
    uint32_t maxTrim,
    uint32_t maxSkip,
    uint32_t minAlignedMarkerCount)
{
    cout << timestamp << "createConflictReadGraph begins." << endl;

    // Check that we have what we need.
    // The code as written only supports the directed read graph.
    SHASTA_ASSERT(directedReadGraph.edges.isOpen);
    SHASTA_ASSERT(directedReadGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(directedReadGraph.edgesByTarget.isOpen());
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the number of threads.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store criteria for the induced alignment so all threads can see them.
    createConflictReadGraphData.inducedAlignmentCriteria.maxOffsetSigma = maxOffsetSigma;
    createConflictReadGraphData.inducedAlignmentCriteria.maxTrim = maxTrim;
    createConflictReadGraphData.inducedAlignmentCriteria.maxSkip = maxSkip;
    createConflictReadGraphData.inducedAlignmentCriteria.minAlignedMarkerCount =
        minAlignedMarkerCount;

    // Initialize the conflict read graph.
    conflictReadGraph.createNew(largeDataName("ConflictReadGraph"), largeDataPageSize);
    conflictReadGraph.createVertices(readCount());

#if 0
    // Compute leftTrim, rightTrim, longestGap for each vertex.
    setupLoadBalancing(conflictReadGraph.vertices.size(), 10);
    runThreads(&Assembler::createConflictReadGraphThreadFunction1, threadCount);
#endif

    // Write a csv file summarizing vertices of the conflict read graph.
    {
        ofstream csv("ConflictReadGraphVertices.csv");
        csv << "VertexId,OrientedReadId\n";
        for(ConflictReadGraph::VertexId v=0; v<conflictReadGraph.vertices.size(); v++) {
            // const ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(v);
            csv <<
                v << "," <<
                ConflictReadGraph::getOrientedReadId(v) <<  "\n";
        }
    }

    // Add edges.
    conflictReadGraph.edges.reserve(10 * readCount());
    setupLoadBalancing(readCount(), 1);
    runThreads(&Assembler::createConflictReadGraphThreadFunction2, threadCount);
    conflictReadGraph.edges.unreserve();
    conflictReadGraph.computeConnectivity();
    conflictReadGraph.writeGraphviz("ConflictReadGraph.dot");


    // Write a csv file summarizing edges of the conflict read graph.
    {
        ofstream csv("ConflictReadGraphEdges.csv");
        csv << "VertexId0,VertexId1,OrientedReadId0,OrientedReadId1,\n";
        for(ConflictReadGraph::EdgeId e=0; e<conflictReadGraph.edges.size(); e++) {
            const ConflictReadGraph::VertexId v0 = conflictReadGraph.v0(e);
            const ConflictReadGraph::VertexId v1 = conflictReadGraph.v1(e);
            csv <<
                v0 << "," <<
                v1 << "," <<
                ConflictReadGraph::getOrientedReadId(v0) << "," <<
                ConflictReadGraph::getOrientedReadId(v1) << "\n";
        }
    }



    cout << "The conflict read graph has " <<
        conflictReadGraph.vertices.size() << " vertices and " <<
        conflictReadGraph.edges.size() << " edges" << endl;
    cout << timestamp << "createConflictReadGraph ends." << endl;
}



void Assembler::accessConflictReadGraph()
{
    conflictReadGraph.accessExistingReadWrite(largeDataName("ConflictReadGraph"));

}


#if 0
void Assembler::createConflictReadGraphThreadFunction1(size_t threadId)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices in this batch.
        for(ConflictReadGraph::VertexId v=begin; v!=end; v++) {
            ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(v);

            // The MarkerId of the first marker for this oriented read.
            const MarkerId firstMarkerId =  markers.begin(v) - markers.begin();

            // The number of markers in this oriented read.
            const uint32_t markerCount = uint32_t(markers.size(v));

            // Compute leftTrim for this vertex.
            bool done = false;
            for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
                const MarkerId markerId = firstMarkerId + ordinal;
                const MarkerGraph::CompressedVertexId compressedVertexId =
                    markerGraph.vertexTable[markerId];
                if(compressedVertexId != MarkerGraph::invalidCompressedVertexId) {
                    vertex.leftTrim = ordinal;
                    done = true;
                    break;
                }
            }
            if(not done) {
                vertex.leftTrim = markerCount;
            }

            // Compute rightTrim for this vertex;
            done = false;
            for(uint32_t ordinal=markerCount-1; ; ordinal--) {
                const MarkerId markerId = firstMarkerId + ordinal;
                const MarkerGraph::CompressedVertexId compressedVertexId =
                    markerGraph.vertexTable[markerId];
                if(compressedVertexId != MarkerGraph::invalidCompressedVertexId) {
                    vertex.rightTrim = markerCount - 1 - ordinal;
                    done = true;
                    break;
                }

                if(ordinal == 0) {
                    break;
                }
            }
            if(not done) {
                vertex.rightTrim = markerCount;
            }

            // Compute longestGap for this vertex;
            vertex.longestGap = markerCount;
            if(vertex.leftTrim + vertex.rightTrim < markerCount-1) {
                vertex.longestGap = 0;
                uint32_t last = vertex.leftTrim;
                for(uint32_t ordinal=vertex.leftTrim+1; ordinal != markerCount-vertex.rightTrim; ordinal++) {
                    const MarkerId markerId = firstMarkerId + ordinal;
                    const MarkerGraph::CompressedVertexId compressedVertexId =
                        markerGraph.vertexTable[markerId];
                    if(compressedVertexId == MarkerGraph::invalidCompressedVertexId) {
                        vertex.longestGap = max(vertex.longestGap, ordinal-last);
                    } else {
                        last = ordinal;
                    }
                }
                if(not(vertex.leftTrim + vertex.rightTrim + vertex.longestGap <= markerCount-2)) {
                    cout <<
                        ConflictReadGraph::getOrientedReadId(v) << " " <<
                        vertex.leftTrim << " " <<
                        vertex.rightTrim << " " <<
                        vertex.longestGap << " " <<
                        markerCount << endl;
                }
                SHASTA_ASSERT(vertex.leftTrim + vertex.rightTrim + vertex.longestGap <= markerCount-2);
            }


            // Set the hasLongGap field.
            vertex.hasLongGap =
                (
                vertex.longestGap >
                createConflictReadGraphData.inducedAlignmentCriteria.maxSkip
                );
        }
    }
}
#endif


void Assembler::createConflictReadGraphThreadFunction2(size_t threadId)
{
    const InducedAlignmentCriteria inducedAlignmentCriteria =
        createConflictReadGraphData.inducedAlignmentCriteria;

    // Work areas for addConflictGraphEdges.
    vector<OrientedReadId> conflictCandidates;
    vector<OrientedReadId> conflictingOrientedReads;
    vector<InducedAlignment> inducedAlignments;
    vector<uint64_t> work;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads in this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); ++readId) {
            // cout << "Working on read " << readId << endl;
            addConflictGraphEdges(
                readId,
                inducedAlignmentCriteria,
                conflictCandidates,
                conflictingOrientedReads,
                inducedAlignments,
                work);
        }
    }

}



// This creates edges of the conflict read graph where
// the lower numbered read is readId0.
// It adds the edges to the conflict read graph directly, under mutex protection.
// This should not create significant contention as adding edges to the
// graph is most of the times much faster than computing them.
void Assembler::addConflictGraphEdges(
    ReadId readId0,
    const InducedAlignmentCriteria& inducedAlignmentCriteria,

    // Work areas.
    vector<OrientedReadId>& conflictCandidates,
    vector<OrientedReadId>& conflictingOrientedReads,
    vector<InducedAlignment>& inducedAlignments,
    vector<uint64_t>& work)
{

    // Put this read on strand 0.
    // When adding edges to the conflict read graph, we will make sure
    // to also add the reverse complemented edge.
    const OrientedReadId orientedReadId0(readId0, 0);
    // const ConflictReadGraph::VertexId v0 = ConflictReadGraph::getVertexId(orientedReadId0);
    // const ConflictReadGraphVertex& vertex0 = conflictReadGraph.getVertex(v0);

#if 0
    // If the vertex corresponding to this read has a long gap don't add any edges.
    if(vertex0.hasLongGap) {
        return;
    }
#endif


    // Find conflict candidates for orientedReadId0.
    // These are OrientedReadId's that share at least one marker graph vertex
    // with orientedReadId0.
    // To do this, we loop over all markers of orientedReadId0.
    conflictCandidates.clear();
    const MarkerId firstMarkerId = markers.begin(orientedReadId0.getValue()) - markers.begin();
    const uint32_t markerCount = uint32_t(markers.size(orientedReadId0.getValue()));
    for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
        const MarkerId markerId0 = firstMarkerId + ordinal;

        // Find the vertex that this marker is on.
        const MarkerGraph::CompressedVertexId compressedVertexId =
            markerGraph.vertexTable[markerId0];

        // If this marker is not on a marker graph vertex, skip.
        if(compressedVertexId == MarkerGraph::invalidCompressedVertexId) {
            continue;
        }

        // Loop over all markers on this vertex.
        const span<MarkerId> vertexMarkers =
            markerGraph.getVertexMarkerIds(compressedVertexId);
        for(const MarkerId markerId1: vertexMarkers) {

            // Skip the marker that we started from.
            if(markerId1 == markerId0) {
                continue;
            }

            // Find the oriented read on this marker.
            OrientedReadId orientedReadId1;
            uint32_t ordinal1;
            tie(orientedReadId1, ordinal1) = findMarkerId(markerId1);

            // Only consider markers with readId1>readId0.
            if(orientedReadId1.getReadId() <= readId0) {
                continue;
            }

            // Add this oriented read to our conflict candidates.
           conflictCandidates.push_back(orientedReadId1);
        }
    }
    deduplicate(conflictCandidates);



    // Remove conflict candidates that correspond to an edge of the read graph.
    // For those we already have a good alignment.
    // Also remove conflict candidates in which the other vertex
    // has a long gap.
    auto itA = conflictCandidates.begin();
    auto itB = itA;
    for(; itA!=conflictCandidates.end(); ++itA) {
        const OrientedReadId orientedReadId1 = *itA;
        const DirectedReadGraph::VertexId v0 = orientedReadId0.getValue();
        const DirectedReadGraph::VertexId v1 = orientedReadId1.getValue();
        const bool forwardExists = directedReadGraph.findEdge(v0, v1)
            != DirectedReadGraph::invalidEdgeId;
        const bool backwardExists = directedReadGraph.findEdge(v1, v0)
            != DirectedReadGraph::invalidEdgeId;
#if 0
        const bool hasLongGap =
            conflictReadGraph.getVertex(v1).hasLongGap;
#endif

        if(/*hasLongGap or */forwardExists or backwardExists) {
            continue;
        } else {
            *itB++ = orientedReadId1;
        }
    }
    conflictCandidates.resize(itB - conflictCandidates.begin());



    // Compute induced alignments between orientedReadId0 and these conflict candidates.
    computeInducedAlignments(
        orientedReadId0,
        conflictCandidates,
        inducedAlignments
    );
    SHASTA_ASSERT(inducedAlignments.size() == conflictCandidates.size());

    // Find which of the induced alignments are bad.
    conflictingOrientedReads.clear();
    // const uint32_t markerCount0 = uint32_t(markers.size(orientedReadId0.getValue()));
    for(uint64_t i=0;i<inducedAlignments.size(); i++) {
        const OrientedReadId orientedReadId1 = conflictCandidates[i];
        // const ConflictReadGraph::VertexId v1 = ConflictReadGraph::getVertexId(orientedReadId1);
        // const ConflictReadGraphVertex& vertex1 = conflictReadGraph.getVertex(v1);
        // const uint32_t markerCount1 = uint32_t(markers.size(orientedReadId1.getValue()));
        if(inducedAlignments[i].indicatesConflict(inducedAlignmentCriteria)) {
            conflictingOrientedReads.push_back(orientedReadId1);
        }
#if 0
        // std::lock_guard<std::mutex> lock(mutex);
        // cout << "Checking induced alignment of " << orientedReadId0 << " " << orientedReadId1 << endl;
        if(not inducedAlignments[i].evaluate(
            markerCount0,
            markerCount1,
            vertex0.leftTrim,
            vertex0.rightTrim,
            vertex1.leftTrim,
            vertex1.rightTrim,
            inducedAlignmentCriteria)) {
            conflictingOrientedReads.push_back(orientedReadId1);
        }
#endif
#if 0
        if(not inducedAlignments[i].evaluate(
            markerCount0,
            markerCount1,
            inducedAlignmentCriteria)) {
            conflictingOrientedReads.push_back(orientedReadId1);
        }
#endif
#if 0
        // This also takes into account the presence or absence of marker graph vertices.
        if(not evaluateInducedAlignment(
            orientedReadId0,
            orientedReadId1,
            inducedAlignments[i],
            inducedAlignmentCriteria,
            work)) {
            conflictingOrientedReads.push_back(orientedReadId1);
        }
#endif
    }
    // cout << "Counts: " << conflictCandidates.size () << " " << conflictingOrientedReads.size() << endl;


    // Add edges to the conflict graph.
    {
        // Find the vertices corresponding to the first read.
        using VertexId = ConflictReadGraph::VertexId;
        const VertexId vertexId0 =
            ConflictReadGraph::getVertexId(orientedReadId0);
        OrientedReadId orientedReadId0ReverseComplement = orientedReadId0;
        orientedReadId0ReverseComplement.flipStrand();
        const VertexId vertexId0ReverseComplement =
            ConflictReadGraph::getVertexId(orientedReadId0ReverseComplement);

        std::lock_guard<std::mutex> lock(mutex);
        for(const OrientedReadId orientedReadId1: conflictingOrientedReads) {

            // Find the vertices corresponding to the second read.
            const VertexId vertexId1 =
                ConflictReadGraph::getVertexId(orientedReadId1);
            OrientedReadId orientedReadId1ReverseComplement = orientedReadId1;
            orientedReadId1ReverseComplement.flipStrand();
            const VertexId vertexId1ReverseComplement =
                ConflictReadGraph::getVertexId(orientedReadId1ReverseComplement);

            // Add the edges.
            conflictReadGraph.addEdge(
                vertexId0,
                vertexId1,
                ConflictReadGraphEdge());
            conflictReadGraph.addEdge(
                vertexId0ReverseComplement,
                vertexId1ReverseComplement,
                ConflictReadGraphEdge());
        }

    }
}



/*******************************************************************************

This removes some vertices from the conflict read graph,
and all edges incident to those vertices.
It also does approximate coloring of the resulting conflict graph.

The vertices and edges are not actually removed,
but their wasRemoved flag is set.
When marking conflict edges in the read graph,
this information is used as follows:
- All read graph edges incident to read graph vertices corresponding to
  conflict read graph vertices that were removed are marked as conflict edges.
  So the oriented reads corresponding to those vertices are
  effectively excluded from the assembly.
- All conflict read graph edges that were removed are ignored.

The goal here is to remove a small number of pathological vertices
with lots of meaningless conflicts.

The current algorithm recursively removes articulation points of
the conflict read graph. This is done by creating a "dynamic"
version of the conflict read graph that uses the Boost library.

*******************************************************************************/

void Assembler::cleanupConflictReadGraph()
{
    using VertexId = ConflictReadGraph::VertexId;
    using EdgeId = ConflictReadGraph::EdgeId;
    using vertex_descriptor = DynamicConflictReadGraph::vertex_descriptor;
    using edge_descriptor = DynamicConflictReadGraph::edge_descriptor;

    // Check that we have what we need.
    SHASTA_ASSERT(conflictReadGraph.isOpen());

    // Create the Boost version of the conflict read graph.
    DynamicConflictReadGraph graph(conflictReadGraph);
    SHASTA_ASSERT(num_vertices(graph) == conflictReadGraph.vertices.size());
    SHASTA_ASSERT(num_edges(graph) == conflictReadGraph.edges.size());
    cout << "The initial conflict graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;



    // Recursively remove articulation points.
    std::map<edge_descriptor, uint64_t> componentMap;
    makeBiconnected(graph, componentMap);
    cout << "After recursive removal of articulation points, the resulting biconnected conflict graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;



    // At this point the graph has no articulation points,
    // which means that the connected components are all biconnected,
    // and therefore the connected components and the biconnected
    // components are the same.
    // We do approximate coloring of each component.
    // Right now the coloring is not used for any purpose.

    // Gather the components.
    std::map<uint64_t, std::set<vertex_descriptor> > components;
    BGL_FORALL_EDGES(e, graph, DynamicConflictReadGraph) {
        const uint64_t componentId = componentMap[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        components[componentId].insert(v0);
        components[componentId].insert(v1);
    }



    // Do approximate coloring of each component.
    uint64_t componentId = 0;
    for(const auto& p: components) {
        const std::set<vertex_descriptor>& componentSet = p.second;
        if(componentSet.size() == 1) {
            continue;
        }
        vector<vertex_descriptor> component;
        copy(componentSet.begin(), componentSet.end(), back_inserter(component));

        // Create a graph representing this component.
        using ComponentGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
        ComponentGraph componentGraph(component.size());
        for(uint64_t iv0=0; iv0<component.size(); iv0++) {
            const vertex_descriptor v0 = component[iv0];
            BGL_FORALL_ADJ(v0, v1, graph, DynamicConflictReadGraph) {
                if(v0 < v1) {
                    const uint64_t iv1 = lower_bound(component.begin(), component.end(), v1) - component.begin();
                    SHASTA_ASSERT(iv1 < component.size());
                    boost::add_edge(iv0, iv1, componentGraph);
                }
            }
        }

        vector<uint64_t> color(component.size());
        approximateColoring(componentGraph, color);
        for(uint64_t i=0; i<component.size(); i++) {
            auto& vertex = graph[component[i]];
            vertex.componentId = componentId;
            vertex.color = color[i];
        }

        vector<uint64_t> colorFrequency;
        deduplicateAndCount(color, colorFrequency);
        cout << "Used " << color.size() << " colors to color a component with " <<
            component.size() << " vertices:" << endl;
        for(uint64_t i=0; i<color.size(); i++) {
            cout << "Color " << color[i] << ": " << colorFrequency[i] << " vertices." << endl;
        }

        componentId++;
    }
    graph.writeGraphviz("ConflictReadGraph-0.dot");



    // Propagate information in the DynamicConflictReadGraph to the ConflictReadGraph.

    // Set the wasRemoved flag for all vertices and edges in the original conflict graph.
    for(VertexId vertexId=0; vertexId<conflictReadGraph.vertices.size(); vertexId++) {
        conflictReadGraph.getVertex(vertexId).wasRemoved = true;
    }
    for(EdgeId edgeId=0; edgeId<conflictReadGraph.edges.size(); edgeId++) {
        conflictReadGraph.getEdge(edgeId).wasRemoved = true;
    }

    // Now clear the wasRemoved flag for the surviving vertices and edges.
    // Also propagate coloring information.
    BGL_FORALL_VERTICES(v, graph, DynamicConflictReadGraph) {
        const VertexId vertexId = graph[v].vertexId;
        auto& vertex = conflictReadGraph.getVertex(vertexId);
        vertex.wasRemoved = false;
        if(graph[v].isColored()) {
            vertex.componentId = uint32_t(graph[v].componentId);
            vertex.color = uint32_t(graph[v].color);
        }
    }
    BGL_FORALL_EDGES(e, graph, DynamicConflictReadGraph) {
        const EdgeId edgeId = graph[e].edgeId;
        conflictReadGraph.getEdge(edgeId).wasRemoved = false;
    }



#if 0
    // Analyze alignments to color more vertices.

    // The vertices that can be colored.
    // (Vertexid, componentId, color).
    vector< tuple<VertexId, uint32_t, uint32_t> > colorableVertices;

    // The ambiguos vertices.
    vector<VertexId> ambiguousVertices;

    // Analyze alignments to color more vertices.
    for(VertexId v0=0; v0<conflictReadGraph.vertices.size(); v0++) {
        ConflictReadGraphVertex cVertex0 = conflictReadGraph.getVertex(v0);

        // If this vertex is already colored, skip it.
        if(cVertex0.hasValidColor()) {
            continue;
        }

        // We will look for alignments where v0 is entirely contained in the other read,
        const uint32_t leftTrim0 = cVertex0.leftTrim;
        const uint32_t rightTrim0 = cVertex0.rightTrim;

        // Set to hold pairs(componentId, color) for containing vertices.
        std::set< pair<uint32_t, uint32_t> > colorInfo;

        // Loop over alignments involving this oriented read.
        // Note this involving using all read graph edges, including the
        // ones not marked as "keep".
        for(EdgeId e01: directedReadGraph.outEdges(v0)) {
            const DirectedReadGraphEdge dEdge01 = directedReadGraph.getEdge(e01);

            // If v0 is not contained, skip.
            const uint32_t maxTrim = 100; // **************** EXPOSE WHEN CODE STABILIZES.
            if(dEdge01.alignmentInfo.leftTrim(0) > leftTrim0 + maxTrim) {
                continue;
            }
            if(dEdge01.alignmentInfo.rightTrim(0) > rightTrim0 + maxTrim) {
                continue;
            }

            const VertexId v1 = directedReadGraph.target(e01);
            const ConflictReadGraphVertex cVertex1 = conflictReadGraph.getVertex(v1);
            if(not cVertex1.hasValidColor()) {
                continue;
            }
            colorInfo.insert(make_pair(cVertex1.componentId, cVertex1.color));

            /*
            cout << ConflictReadGraph::getOrientedReadId(v0) <<
                " is contained in " << ConflictReadGraph::getOrientedReadId(v1) <<
                " with componentId " << cVertex1.componentId <<
                " color " << cVertex1.color << endl;
            */
        }
        for(EdgeId e10: directedReadGraph.inEdges(v0)) {
            const DirectedReadGraphEdge dEdge10 = directedReadGraph.getEdge(e10);

            // If v0 is not contained, skip.
            const uint32_t maxTrim = 100; // **************** EXPOSE WHEN CODE STABILIZES.
            if(dEdge10.alignmentInfo.leftTrim(1) > leftTrim0 + maxTrim) {
                continue;
            }
            if(dEdge10.alignmentInfo.rightTrim(1) > rightTrim0 + maxTrim) {
                continue;
            }

            const VertexId v1 = directedReadGraph.source(e10);
            const ConflictReadGraphVertex cVertex1 = conflictReadGraph.getVertex(v1);
            if(not cVertex1.hasValidColor()) {
                continue;
            }
            colorInfo.insert(make_pair(cVertex1.componentId, cVertex1.color));

            /*
            cout << ConflictReadGraph::getOrientedReadId(v0) <<
                " is contained in " << ConflictReadGraph::getOrientedReadId(v1) <<
                " with componentId " << cVertex1.componentId <<
                " color " << cVertex1.color << endl;
            */
        }

        if(colorInfo.size() == 1) {
            // cout << "Unambiguous coloring possible for " << ConflictReadGraph::getOrientedReadId(v0) << endl;
            const auto& p = *colorInfo.begin();
            colorableVertices.push_back(make_tuple(v0, p.first, p.second));
        } else if(colorInfo.size() > 1) {
            ambiguousVertices.push_back(v0);
            // cout << "Ambiguous containments for " << ConflictReadGraph::getOrientedReadId(v0) << endl;
        }
    }
    cout << "After containment analysis, found " <<
        colorableVertices.size() << " vertices that can be colored and " <<
        ambiguousVertices.size() << " ambiguous vertices." << endl;

    // Color them.
    for(const auto& t: colorableVertices) {
        ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(get<0>(t));
        vertex.componentId = get<1>(t);
        vertex.color = get<2>(t);
    }

    // Mark the ambiguous vertices as removed.
    for(const VertexId v: ambiguousVertices) {
        conflictReadGraph.getVertex(v).wasRemoved = true;
        for(const EdgeId e: conflictReadGraph.incidentEdges(v)) {
            conflictReadGraph.getEdge(e).wasRemoved = true;
        }
    }
#endif
}



#if 0
// This colors the ConflictReadGraph by walking the DirectedReadGraph.
void Assembler::colorConflictReadGraph()
{
    const bool debug = false;

    // Check that we have what we need.
    SHASTA_ASSERT(directedReadGraph.isOpen());
    SHASTA_ASSERT(conflictReadGraph.isOpen());

    // Types for vertices and edges of the two graphs we will use.
    using VertexId = DirectedReadGraph::VertexId;
    using EdgeId   = DirectedReadGraph::EdgeId;

    // We are assuming the two graphs use the same VertexId and EdgeId.
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");

    // Initialize all vertex clusterId's to invalid.
    const auto invalid = ConflictReadGraphVertex::invalid;
    const VertexId n = conflictReadGraph.vertices.size();
    for(VertexId vertexId=0; vertexId<n; vertexId++) {
        auto& vertex = conflictReadGraph.getVertex(vertexId);
        vertex.componentId = invalid;
        vertex.color = invalid;
    }


    // Create a table of vertices excluding vertices with long gap and containing:
    // - VertexId.
    // - Number of conflicts (degree in the conflict read graph).
    // - Degree in the read graph, counting only edges flagged as "keep".
    // Sort it by increasing number of conflicts, then by decreasing degree.
    vector<ColorConflictReadGraphData> vertexTable;
    for(VertexId v=0; v<n; v++) {
        if(not conflictReadGraph.getVertex(v).hasLongGap) {
            vertexTable.push_back(
                ColorConflictReadGraphData(v, directedReadGraph, conflictReadGraph));
        }
    }
    sort(vertexTable.begin(), vertexTable.end());



    // Data structures to keep track, at each iteration, of the vertices
    // we already encountered.
    vector<VertexId> encounteredVertices;
    vector<bool> wasEncountered(directedReadGraph.vertices.size(), false);

    // Data structures to keep track, at each iteration, of the vertices
    // that conflict with vertices we already encountered.
    vector<VertexId> forbiddenVertices;
    vector<bool> isForbidden(directedReadGraph.vertices.size(), false);

    // Other data structures used below.
    vector<VertexId> adjacentVertices;
    vector< pair<VertexId, uint64_t > > adjacentVerticesSortedByConflictCount;



    // Iterate over possible starting vertices in the order in which
    // they appear in the vertex table.
    uint64_t vertexTableIndex = 0;
    uint32_t iteration=0;
    for(; ; iteration++) {
        if(vertexTableIndex == vertexTable.size()) {
            break;
        }

        // Find the next vertex that has not yet been colored (assigned to a cluster).
        VertexId startVertexId = vertexTable[vertexTableIndex++].vertexId;
        bool done = false;
        while(true) {
            if(not conflictReadGraph.getVertex(startVertexId).hasValidColor()) {
                break;
            }
            if(vertexTableIndex == n) {
                done = true;
                break;
            }
            startVertexId = vertexTable[vertexTableIndex++].vertexId;
        }
        if(done) {
            break;
        }
        if(debug) {
            cout << "Start iteration " << iteration <<
                " from " << ConflictReadGraph::getOrientedReadId(startVertexId);
            cout << " with number of conflicts " <<
                vertexTable[vertexTableIndex-1].conflictReadGraphDegree <<
                " and kept degree  " <<
                vertexTable[vertexTableIndex-1].directedReadGraphKeptDegree << endl;
        }
        SHASTA_ASSERT(not conflictReadGraph.getVertex(startVertexId).hasValidColor());


        // We use a process similar to a BFS starting at this vertex:
        // - When encountering a vertex that conflicts with another vertex
        //   we already encountered at this iteration, we skip it.
        // - When encountering a vertex that was already colored at a previous
        //   iteration (not just one colored at the current iteration), we skip it.
        // - When we enqueue neighbors of a vertex, we enqueue them in order
        //   of increasing number of conflicting vertices.


        // Initialize the BFS.
        std::priority_queue<ColorConflictReadGraphData,
            vector<ColorConflictReadGraphData>,
            std::greater<ColorConflictReadGraphData> > q;
        q.push(ColorConflictReadGraphData(startVertexId, directedReadGraph, conflictReadGraph));
        wasEncountered[startVertexId] = true;
        encounteredVertices.push_back(startVertexId);
        uint64_t coloredCount = 0;

        // Mark as forbidden the vertices that conflict with startVertex.
        for(const EdgeId e: conflictReadGraph.edgesByVertex[startVertexId]) {
            const VertexId v = conflictReadGraph.otherVertex(e, startVertexId);
            if(not isForbidden[v]) {
                isForbidden[v] = true;
                forbiddenVertices.push_back(v);
                if(debug) {
                    cout << "Marked as forbidden " <<
                        conflictReadGraph.getOrientedReadId(v) << endl;
                }
            }
        }

        // BFS loop.
        while(not q.empty()) {

            // Dequeue a vertex.
            const VertexId v0 = q.top().vertexId;
            const OrientedReadId orientedReadId0 = ConflictReadGraph::getOrientedReadId(v0);
            if(debug) {
                cout << "Queue size " << q.size() << ", dequeued " << orientedReadId0 << " " <<
                    q.top().conflictReadGraphDegree << " " <<
                    q.top().directedReadGraphKeptDegree << endl;
            }
            q.pop();

            // If v0 is now forbidden, skip it.
            if(isForbidden[v0]) {
                // v0 was not forbidden when we enqueud it, but it is forbidden now.
                if(debug) {
                    cout << orientedReadId0 << " skipped because it is now forbidden." << endl;
                }
                continue;
            }

            // Give it a color equal to this iteration.
            SHASTA_ASSERT(not conflictReadGraph.getVertex(v0).hasValidColor());
            conflictReadGraph.getVertex(v0).color = iteration;
            coloredCount++;
            if(debug) {
                cout<< orientedReadId0 << " being assigned to cluster id " << iteration << endl;
            }

            // Gather adjacent vertices.
            directedReadGraph.findKeptAdjacent(v0, adjacentVertices);

            // Loop over adjacent vertices, in this order of increasing number of conflicts.
            // Skip vertices with a long gap.
            for(const VertexId v1: adjacentVertices) {

                // See if this vertex should be skipped.
                if(isForbidden[v1]) {
                    if(debug) {
                        cout << ConflictReadGraph::getOrientedReadId(v1) << " forbidden" << endl;
                    }
                    continue;
                }
                if(conflictReadGraph.vertices[v1].color != invalid) {
                    if(debug) {
                        cout << ConflictReadGraph::getOrientedReadId(v1) << " already assigned to a cluster" << endl;
                    }
                    continue;
                }
                if(wasEncountered[v1]) {
                    if(debug) {
                        cout << ConflictReadGraph::getOrientedReadId(v1) << " already encountered" << endl;
                    }
                    continue;
                }
                if(conflictReadGraph.vertices[v1].hasLongGap) {
                    if(debug) {
                        cout << ConflictReadGraph::getOrientedReadId(v1) << " has long gap" << endl;
                    }
                    continue;
                }

                // We know this vertex is not forbidden and was not already
                // colored at this or the previous iteration, so we can enqueue it now.
                SHASTA_ASSERT(not wasEncountered[v1]);
                wasEncountered[v1] = true;
                encounteredVertices.push_back(v1);
                const ColorConflictReadGraphData data1(v1, directedReadGraph, conflictReadGraph);
                q.push(data1);
                if(debug) {
                    cout << "Enqueued " << conflictReadGraph.getOrientedReadId(v1) << " " <<
                        data1.conflictReadGraphDegree << " " <<
                        data1.directedReadGraphKeptDegree << endl;
                }

                // Mark as forbidden the vertices that conflict with v1.
                for(const EdgeId e12: conflictReadGraph.edgesByVertex[v1]) {
                    const VertexId v2 = conflictReadGraph.otherVertex(e12, v1);
                    if(not isForbidden[v2]) {
                        isForbidden[v2] = true;
                        forbiddenVertices.push_back(v2);
                        if(debug) {
                            cout << "Marked as forbidden " <<
                                conflictReadGraph.getOrientedReadId(v2) << endl;
                        }
                    }
                }
            }
        }



        // Clean up data structures to prepare them for the next iteration.
        for(const VertexId v: forbiddenVertices) {
            SHASTA_ASSERT(isForbidden[v]);
            isForbidden[v] = false;
        }
        forbiddenVertices.clear();
        for(const VertexId v: encounteredVertices) {
            SHASTA_ASSERT(wasEncountered[v]);
            wasEncountered[v] = false;
        }
        encounteredVertices.clear();

        if(debug) {
            cout << "Iteration " << iteration << " colored " <<
                coloredCount << " vertices." << endl;
        }
    }



    // Check that all vertices were assigned to clusters except for those
    // that have a long gap.
    uint32_t longGapCount = 0;
    for(VertexId v=0; v<n; v++) {
        const ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(v);
        if(not vertex.hasLongGap) {
            SHASTA_ASSERT(conflictReadGraph.getVertex(v).hasValidClusterId());
        } else {
            ++longGapCount;
        }
    }

    // Store the vertices in each cluster.
    vector< vector<VertexId> > clusters(iteration);
    for(VertexId v=0; v<n; v++) {
        const ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(v);
        if(not vertex.hasLongGap) {
            clusters[conflictReadGraph.getVertex(v).clusterId].push_back(v);
        }
    }




    // Uncolor small clusters, and write a summary of the clusters we found.
    const uint64_t clusterSizeThreshold = 10;  // ********************* EXPOSE WHEN CODE STABILIZES.
    uint64_t smallClusterVertexCount = 0;
    for(uint64_t clusterId=0; clusterId<clusters.size(); clusterId++) {
        const vector<VertexId>& cluster = clusters[clusterId];
        if(cluster.size() >= clusterSizeThreshold) {
            cout << "Cluster " << clusterId << " has " << cluster.size() << " vertices." << endl;
        } else {
            for(const VertexId v: cluster) {
                conflictReadGraph.getVertex(v).clusterId = invalid;
                ++smallClusterVertexCount;
            }
        }
    }
    cout << "Number of vertices in small, discarded clusters: "
        << smallClusterVertexCount << endl;
    cout << "Number of vertices with a long gap, not assigned to any cluster: "
        << longGapCount << endl;
    cout << "Total number of vertices " << conflictReadGraph.vertices.size() << endl;



    // Write out a statistics file with a line for each edge
    // in the conflict read graph.
    ofstream csv("ColoringStatistics.csv");
    csv << "OrientedReadId0,OrientedREadId1,ClusterId0,ClusterId1,MinClusterId,MaxClusterId,ClusterSize0,ClusterSize1\n";
    for(EdgeId e=0; e<conflictReadGraph.edges.size(); e++) {
        const VertexId v0 = conflictReadGraph.v0(e);
        const VertexId v1 = conflictReadGraph.v1(e);
        const uint64_t clusterId0 = conflictReadGraph.getVertex(v0).clusterId;
        const uint64_t clusterId1 = conflictReadGraph.getVertex(v1).clusterId;
        /*
        if(clusterId0 == clusterId1) {
            cout << "Inconsistent coloring " <<
                ConflictReadGraph::getOrientedReadId(v0) << " " <<
                ConflictReadGraph::getOrientedReadId(v1) << endl;
        }
        */
        csv <<
            ConflictReadGraph::getOrientedReadId(v0) << "," <<
            ConflictReadGraph::getOrientedReadId(v1) << "," <<
            clusterId0 << "," <<
            clusterId1 << "," <<
            min(clusterId0, clusterId1) << "," <<
            max(clusterId0, clusterId1) << ",";
        if(clusterId0!=invalid and clusterId1!=invalid) {
            SHASTA_ASSERT(clusterId0 != clusterId1);
            csv <<
                clusters[clusterId0].size() << "," <<
                clusters[clusterId1].size();
        }
        csv << "\n";
    }



    // Create a map that contains all conflict graph edges between each pair of clusters.
    std::map< pair<uint64_t, uint64_t>, vector<EdgeId> > conflictGraphEdgesBetweenClusters;
    for(EdgeId e=0; e<conflictReadGraph.edges.size(); e++) {
        const VertexId v0 = conflictReadGraph.v0(e);
        const VertexId v1 = conflictReadGraph.v1(e);
        uint64_t clusterId0 = conflictReadGraph.getVertex(v0).clusterId;
        uint64_t clusterId1 = conflictReadGraph.getVertex(v1).clusterId;
        if(clusterId0==invalid || clusterId1==invalid) {
            continue;
        }
        SHASTA_ASSERT(clusterId0 != clusterId1);
        if(clusterId1 < clusterId0) {
            swap(clusterId0, clusterId1);
        }
        conflictGraphEdgesBetweenClusters[make_pair(clusterId0,clusterId1)].push_back(e);
    }



    // Check pairs of clusters with conflicts between them.
    for(const auto& p: conflictGraphEdgesBetweenClusters) {

        const uint64_t clusterId0 = p.first.first;
        const vector<VertexId>& cluster0 = clusters[clusterId0];
        const uint64_t clusterSize0 = cluster0.size();

        const uint64_t clusterId1 = p.first.second;
        const vector<VertexId>& cluster1 = clusters[clusterId1];
        const uint64_t clusterSize1 = cluster1.size();

        const vector<EdgeId>& conflictEdges =p.second;
        cout << conflictEdges.size() << " conflicts between clusters " <<
            clusterId0 << " " << clusterId1 <<
            " of sizes " << clusterSize0 << " " << clusterSize1 << endl;
    }

}
#endif



// This uses the coloring of the conflict graph.
// It walks the read graph while avoiding conflicts.
void Assembler::markDirectedReadGraphConflictEdges1()
{
    // Check that we have what we need.
    SHASTA_ASSERT(directedReadGraph.isOpen());
    SHASTA_ASSERT(conflictReadGraph.isOpen());

    // Loop over all edges of the directed read graph.
    uint64_t invalidEdgeCount = 0;
    uint64_t invalidKeptEdgeCount = 0;
    uint64_t keptEdgeCount = 0;
    for(DirectedReadGraph::EdgeId edgeId=0; edgeId<directedReadGraph.edges.size(); edgeId++) {
        DirectedReadGraphEdge& edge = directedReadGraph.getEdge(edgeId);

        // Get the vertices of the DirectedReadGraph..
        const DirectedReadGraph::VertexId v0 = directedReadGraph.source(edgeId);
        const DirectedReadGraph::VertexId v1 = directedReadGraph.target(edgeId);

        // Get the corresponding OrientedReadId's.
        const OrientedReadId orientedReadId0 = OrientedReadId(OrientedReadId::Int(v0));
        const OrientedReadId orientedReadId1 = OrientedReadId(OrientedReadId::Int(v1));

        // Get the corresponding vertices of the ConflictReadGraph.
        const ConflictReadGraph::VertexId u0 = ConflictReadGraph::getVertexId(orientedReadId0);
        const ConflictReadGraph::VertexId u1 = ConflictReadGraph::getVertexId(orientedReadId1);
        const ConflictReadGraphVertex& cVertex0 = conflictReadGraph.getVertex(u0);
        const ConflictReadGraphVertex& cVertex1 = conflictReadGraph.getVertex(u1);

        // With current numbering, the vertex ids should be the same.
        SHASTA_ASSERT(u0 == v0);
        SHASTA_ASSERT(u1 == v1);

        // Figure out if this a conflict edge.
        edge.isConflict =
            (not cVertex0.hasValidColor()) or
            (not cVertex1.hasValidColor()) or
            (cVertex0.color != cVertex1.color);

        if(edge.isConflict) {
            ++invalidEdgeCount;
        }

        if(edge.keep) {
            ++keptEdgeCount;
            if(edge.isConflict) {
                ++invalidKeptEdgeCount;
            }
        }
    }

    cout << "Directed read graph edge counts:" << endl;
    cout << "    Total " << directedReadGraph.edges.size() << endl;
    cout << "    Kept for marker graph creation " << keptEdgeCount << endl;
    cout << "    Marked as conflict " << invalidEdgeCount << endl;
    cout << "    Kept for marker graph creation and marked as conflict " << invalidKeptEdgeCount << endl;
}



Assembler::ColorConflictReadGraphData::ColorConflictReadGraphData(
    DirectedReadGraph::VertexId vertexId,
    const DirectedReadGraph& directedReadGraph,
    const ConflictReadGraph& conflictReadGraph) :
    vertexId(vertexId),
    conflictReadGraphDegree(conflictReadGraph.degree(vertexId)),
    directedReadGraphKeptDegree(directedReadGraph.keptDegree(vertexId))
{
}



// This does not use the coloring of conflict graph.
// If gradually adds edges to the read graph, making sure
// to not create conflicts in read graph neighborhoods
// of the specified radius.
void Assembler::markDirectedReadGraphConflictEdges2(int radius)
{
    const bool debug = false;

    // Check that we have what we need.
    SHASTA_ASSERT(directedReadGraph.isOpen());
    SHASTA_ASSERT(conflictReadGraph.isOpen());

    // Types for vertices and edges of the two graphs we will use.
    using VertexId = DirectedReadGraph::VertexId;
    using EdgeId   = DirectedReadGraph::EdgeId;

    // We are assuming the two graphs use the same VertexId and EdgeId.
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");

    // Begin by marking all read graph edges as "isConflict".
    for(EdgeId e=0; e<directedReadGraph.edges.size(); e++) {
        DirectedReadGraphEdge& edge = directedReadGraph.getEdge(e);
        edge.isConflict = 1;
    }


    // Gather all read graph edges that:
    // 1. Are marked "keep"
    // AND
    // 2. Do not involve vertices that were removed from the ConflictReadGraph.
    // For each edge, also store the number aligned markers,
    // then sort by decreasing number of aligned markers.
    vector< pair<VertexId, uint32_t> > edges;
    for(EdgeId e=0; e<directedReadGraph.edges.size(); e++) {
        const DirectedReadGraphEdge& edge = directedReadGraph.getEdge(e);

        // If the edge is not marked "keep", skip it.
        if(not edge.keep) {
            continue;
        }

        // Check that the vertices are were not removed from the ConflictReadGraph.
        const VertexId v0 = directedReadGraph.source(e);
        if(conflictReadGraph.getVertex(v0).wasRemoved) {
            continue;
        }
        const VertexId v1 = directedReadGraph.target(e);
        if(conflictReadGraph.getVertex(v1).wasRemoved) {
            continue;
        }

        edges.push_back(make_pair(e, edge.alignmentInfo.markerCount));
    }
    sort(edges.begin(), edges.end(),
        OrderPairsBySecondOnlyGreater<VertexId, uint32_t>());



    // Create an edge filter that allows edges
    // marked as "keep" and not marked as "isConflict".
    const DirectedReadGraph::EdgeFilter edgeFilter(
        0,
        std::numeric_limits<uint64_t>::max(),
        0.,
        false,
        true);



    // Loop over these edges in this order. Only add an edge to read
    // graph (by clearing the "isConflict" flag) if adding it
    // does not introduce a conflict within the specified radius of the edge.
    std::map<VertexId, uint64_t> neighborMap;
    vector<VertexId> neighborhood;
    for(const auto& p: edges) {
        const EdgeId e = p.first;
        SHASTA_ASSERT(directedReadGraph.getEdge(e).isConflict == 1);
        const VertexId v0 = directedReadGraph.source(e);
        const VertexId v1 = directedReadGraph.target(e);
        if(debug) {
            cout << "Checking for conflicts read graph edge " <<
                ConflictReadGraph::getOrientedReadId(v0) << "->" <<
                ConflictReadGraph::getOrientedReadId(v1) << endl;
        }

        neighborhood.clear();
        directedReadGraph.findNeighborhood(v0, radius, edgeFilter, true, true, 0., neighborMap);
        for(const auto& p: neighborMap) {
            neighborhood.push_back(p.first);
        }
        directedReadGraph.findNeighborhood(v1, radius, edgeFilter, true, true, 0., neighborMap);
        for(const auto& p: neighborMap) {
            neighborhood.push_back(p.first);
        }

        // Deduplicate and sort the vertices in this neighborhood.
        deduplicate(neighborhood);

        // Look for edges in the conflict read graph involving vertices in this neighborhood.
        bool conflictEdgeWasFound = false;
        for(const VertexId u0: neighborhood) {
            const span<EdgeId> edges0 = conflictReadGraph.incidentEdges(u0);
            for(const EdgeId e01: edges0) {
                const VertexId u1 = conflictReadGraph.otherVertex(e01, u0);
                if(binary_search(neighborhood.begin(), neighborhood.end(), u1)) {
                    conflictEdgeWasFound = true;

                    if(debug) {
                        cout << "Read graph edge " <<
                            ConflictReadGraph::getOrientedReadId(v0) << "->" <<
                            ConflictReadGraph::getOrientedReadId(v1) <<
                            " marked as conflict because of conflict graph edge " <<
                            ConflictReadGraph::getOrientedReadId(u0) << "->" <<
                            ConflictReadGraph::getOrientedReadId(u1) << endl;
                    }

                    break;
                }
            }
        }

        // If we did not find a conflict, we can add this edge to the read
        // graph by clearing its "isConflict" flag.
        if(not conflictEdgeWasFound) {
            directedReadGraph.getEdge(e).isConflict = 0;
        }
    }

}



void Assembler::markDirectedReadGraphConflictEdges3(int radius)
{
    const bool debug = false;

    // Check that we have what we need.
    SHASTA_ASSERT(directedReadGraph.isOpen());
    SHASTA_ASSERT(conflictReadGraph.isOpen());

    // Types for vertices and edges of the two graphs we will use.
    using VertexId = DirectedReadGraph::VertexId;
    using EdgeId   = DirectedReadGraph::EdgeId;

    // We are assuming the two graphs use the same VertexId and EdgeId.
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");

    // Begin by marking all read graph edges as "isConflict".
    for(EdgeId e=0; e<directedReadGraph.edges.size(); e++) {
        DirectedReadGraphEdge& edge = directedReadGraph.getEdge(e);
        edge.isConflict = 1;
    }


    // Gather all read graph edges that:
    // 1. Are marked "keep"
    // AND
    // 2. Do not involve vertices that were removed from the ConflictReadGraph.
    // For each edge, also store the number aligned markers,
    // then sort by decreasing number of aligned markers.
    vector< pair<VertexId, uint32_t> > edges;
    for(EdgeId e=0; e<directedReadGraph.edges.size(); e++) {
        const DirectedReadGraphEdge& edge = directedReadGraph.getEdge(e);

        // If the edge is not marked "keep", skip it.
        if(not edge.keep) {
            continue;
        }

        // Check that the vertices are were not removed from the ConflictReadGraph.
        const VertexId v0 = directedReadGraph.source(e);
        if(conflictReadGraph.getVertex(v0).wasRemoved) {
            continue;
        }
        const VertexId v1 = directedReadGraph.target(e);
        if(conflictReadGraph.getVertex(v1).wasRemoved) {
            continue;
        }

        edges.push_back(make_pair(e, edge.alignmentInfo.markerCount));
    }
    sort(edges.begin(), edges.end(),
        OrderPairsBySecondOnlyGreater<VertexId, uint32_t>());



    // Create an edge filter that allows edges
    // marked as "keep" and not marked as "isConflict".
    const DirectedReadGraph::EdgeFilter edgeFilter(
        0,
        std::numeric_limits<uint64_t>::max(),
        0.,
        false,
        true);



    // Loop over these edges in this order. Only add an edge to the read
    // graph (by clearing the "isConflict" flag) if adding it
    // does not introduce a path, within the specified radius of the edge,
    // between vertices belonging to the same connected component and different colors
    // of the conflict read graph.
    std::map<VertexId, uint64_t> neighborMap;
    vector<VertexId> neighbors;
    for(const auto& p: edges) {
        const EdgeId e = p.first;
        SHASTA_ASSERT(directedReadGraph.getEdge(e).isConflict == 1);
        const VertexId v0 = directedReadGraph.source(e);
        const VertexId v1 = directedReadGraph.target(e);
        if(debug) {
            cout << "Checking for conflicts read graph edge " <<
                ConflictReadGraph::getOrientedReadId(v0) << "->" <<
                ConflictReadGraph::getOrientedReadId(v1) << endl;
        }

        // Compute the neighbors of v0 and v1, within the specified radius.
        neighbors.clear();
        directedReadGraph.findNeighborhood(v0, radius, edgeFilter, true, true, 0., neighborMap);
        for(const auto& p: neighborMap) {
            neighbors.push_back(p.first);
        }
        directedReadGraph.findNeighborhood(v1, radius, edgeFilter, true, true, 0., neighborMap);
        for(const auto& p: neighborMap) {
            neighbors.push_back(p.first);
        }
        deduplicate(neighbors);

        if(debug) {
            cout << "Neighbors: ";
            for(const VertexId vertexId: neighbors) {
                cout << ConflictReadGraph::getOrientedReadId(vertexId);
                const ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(vertexId);
                SHASTA_ASSERT(!vertex.wasRemoved);
                if(vertex.hasValidColor()) {
                    cout << vertex.componentId << " " << vertex.color << endl;
                }
            }
        }

        // Find the colors present for each component.
        // The colorTable is indexed by component id and contains
        // the set of colors present in the neighbors for each component.
        std::map<uint32_t, std::set<uint32_t> > colorTable;
        for(const VertexId vertexId: neighbors) {
            const ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(vertexId);
            SHASTA_ASSERT(!vertex.wasRemoved);
            if(vertex.hasValidColor()) {
                colorTable[vertex.componentId].insert(vertex.color);
            }
        }

        // Check for conflicts.
        // A conflict occurs if any component has more than one color present.
        bool conflictWasFound = false;
        for(const auto& p: colorTable) {
            if(p.second.size() > 1) {
                conflictWasFound = true;
                break;
            }
        }

        // If we did not find a conflict, we can add this edge to the read
        // graph by clearing its "isConflict" flag.
        if(not conflictWasFound) {
            if(debug)  {
                cout << "No conflict found." << endl;
            }
            directedReadGraph.getEdge(e).isConflict = 0;
        }
    }



    // If an edge is marked as conflict, also mark the reverse complement as conflict.
    for(EdgeId edgeId=0; edgeId<directedReadGraph.edges.size(); edgeId++) {
        const DirectedReadGraphEdge& edge = directedReadGraph.getEdge(edgeId);
        if(edge.isConflict == 1) {
            const EdgeId reverseComplementEdgeId = edge.reverseComplementedEdgeId;
            DirectedReadGraphEdge& reverseComplementEdge =
                directedReadGraph.getEdge(reverseComplementEdgeId);
            reverseComplementEdge.isConflict = 1;
        }
    }

}
