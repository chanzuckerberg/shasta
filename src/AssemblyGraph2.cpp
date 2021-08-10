#include "AssemblyGraph2.hpp"
#include "AssembledSegment.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "computeLayout.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
#include "orderPairs.hpp"
#include "writeGraph.hpp"
using namespace shasta;

#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "fstream.hpp"
#include <limits>
#include <map>
#include <numeric>


// The constructor creates an edge for each linear path
// in the marker graph. Therefore, immediately after construction,
// each edge has a single MarkerGraphPath (no bubbles).
AssemblyGraph2::AssemblyGraph2(
    uint64_t k, // Marker length
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    k(k),
    markers(markers),
    markerGraph(markerGraph)
{

    // Because of the way we write the GFA file (without overlaps),
    // k is required to be even.
    SHASTA_ASSERT((k % 2) == 0);



    // Parameters that should be exposed as options when the
    // code stabilizes.

    // Maximum period to remove bubbles caused by copy number differences.
    const uint64_t maxPeriod = 4;

    // Minimum number of supporting reads for a bubble graph edge
    // to be kept.
    const uint64_t minReadCount = 3;

    // Threshold used to remove bad bubbles.
    const double discordantRatioThreshold = 0.2;

    // Ambiguity threshold for bubble graph edges.
    double ambiguityThreshold = 0.5;

    // Create the assembly graph.
    cout << timestamp << "AssemblyGraph2::create begins." << endl;
    create();

    // Gather parallel edges into bubbles.
    cout << timestamp << "AssemblyGraph2::gatherBubbles begins." << endl;
    gatherBubbles();

    // Store the reads supporting each branch of each edges.
    cout << timestamp << "AssemblyGraph2::storeReadInformation begins." << endl;
    storeReadInformation();

    // Remove bubbles caused by secondary edges.
    cout << timestamp << "AssemblyGraph2::removeSecondaryBubbles begins." << endl;
    removeSecondaryBubbles();

    // Merge adjacent non-bubbles created by the removal of secondary bubbles.
    cout << timestamp << "AssemblyGraph2::merge begins." << endl;
    merge();

    // Assemble sequence.
    cout << timestamp <<"AssemblyGraph2::assemble begins." << endl;
    assemble();

    // Store the sequence to be written to gfa.
    // This transfers common sequence out of bubbles and into adjacent
    // non-bubbles.
    cout << timestamp << "AssemblyGraph2::storeGfaSequence begins." << endl;
    storeGfaSequence();

    // Find bubbles caused by copy number changes in repeats
    // with period up to maxPeriod.
    cout << timestamp << "AssemblyGraph2::findCopyNumberBubbles begins." << endl;
    findCopyNumberBubbles(maxPeriod);

    // Create the bubble graph.
    cout << timestamp << "AssemblyGraph2::createBubbleGraph begins." << endl;
    createBubbleGraph(markers.size()/2);
    cout << "The initial bubble graph has " << num_vertices(bubbleGraph) <<
        " vertices and " << num_edges(bubbleGraph) << " edges." << endl;
    if(false) {
        bubbleGraph.writeGraphviz("BubbleGraph-0.dot");
        bubbleGraph.writeEdgesCsv("BubbleGraphEdges-0.csv");
    }

    // Cleanup the bubble graph.
    // This marks as bad the bubbles corresponding to bubble graph vertices
    // that are removed.
    cout << timestamp << "AssemblyGraph2::cleanupBubbleGraph begins." << endl;
    cleanupBubbleGraph(minReadCount, discordantRatioThreshold, ambiguityThreshold);

    // Compute connected components of the bubble graph.
    bubbleGraph.computeConnectedComponents();

    // Use each connected component of the bubble graph to phase the bubbles.
    cout << timestamp << "AssemblyGraph2::phase begins." << endl;
    phase();

    // Write out what we have.
    cout << timestamp << "AssemblyGraph2 writing GFA output." << endl;
    const bool writeSequence = true;
    writeGfaBothStrands("Assembly-BothStrands", writeSequence);

}



void AssemblyGraph2::cleanupBubbleGraph(
    uint64_t minReadCount,
    double discordantRatioThreshold,
    double ambiguityThreshold)
{
    G& g = *this;
    const bool debug = false;

    // Remove bubble graph edges with low read support.
    bubbleGraph.removeWeakEdges(minReadCount);
    cout << "After removing edges supported by less than " << minReadCount <<
        " reads, the bubble graph has " << num_vertices(bubbleGraph) <<
        " vertices and " << num_edges(bubbleGraph) << " edges." << endl;
    if(debug) {
        bubbleGraph.writeGraphviz("BubbleGraph-1.dot");
        bubbleGraph.writeEdgesCsv("BubbleGraphEdges-1.csv");
    }

    // Remove weak vertices of the bubble graph
    // and flag the corresponding bubbles as bad.
    for(uint64_t iteration=0; ; iteration++) {
        vector<AssemblyGraph2::edge_descriptor> badBubbles;
        bubbleGraph.removeWeakVertices(discordantRatioThreshold, badBubbles);
        if(badBubbles.empty()) {
            break;
        }
        for(const AssemblyGraph2::edge_descriptor e: badBubbles) {
            g[e].isBad = true;
        }
        cout << "After removing " << badBubbles.size() <<
            " weak vertices, the bubble graph has " << num_vertices(bubbleGraph) <<
            " vertices and " << num_edges(bubbleGraph) << " edges." << endl;
        if(debug) {
            bubbleGraph.writeGraphviz("BubbleGraph-Iteration-" + to_string(iteration) + ".dot");
            bubbleGraph.writeEdgesCsv("BubbleGraphEdges-Iteration-" + to_string(iteration) + ".csv");
        }
    }


#if 0
    // Finally, remove edges with high ambiguity.
    vector<BubbleGraph::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        if(bubbleGraph[e].ambiguity() > ambiguityThreshold) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const BubbleGraph::edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, bubbleGraph);
    }
    cout << "Removed " << edgesToBeRemoved.size() <<
        " bubble graph edges with ambiguity greater than " <<
        ambiguityThreshold << endl;
#endif

    if(debug) {
        bubbleGraph.writeGraphviz("BubbleGraph-Final.dot");
        bubbleGraph.writeEdgesCsv("BubbleGraphEdges-Final.csv");
    }
}



// Initial creation of vertices and edges.
void AssemblyGraph2::create()
{

    const bool debug = false;
    ofstream debugOut;
    if(debug) {
        debugOut.open("AssemblyGraph2-Constructor.txt");
    }

    const MarkerGraph::EdgeId edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;



    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(MarkerGraph::EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {
        if(debug) {
            const MarkerGraph::Edge& startEdge = markerGraph.edges[startEdgeId];
            debugOut << "Starting a new path at edge " << startEdgeId << " " <<
                startEdge.source << "->" << startEdge.target << endl;
        }

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId v1 = edge.target;
            const auto outEdges = markerGraph.edgesBySource[v1];
            if(outEdges.size() != 1) {
                break;
            }
            const auto inEdges = markerGraph.edgesByTarget[v1];
            if(inEdges.size() != 1) {
                break;
            }
            edgeId = outEdges[0];
            if(debug) {
                debugOut << "Forward " << edgeId << " " <<
                    edge.source << "->" << edge.target << endl;
            }
            if(edgeId == startEdgeId) {
                isCircular = true;
                if(debug) {
                    cout << "Found a circular edge." << endl;
                    debugOut << "Found a circular edge." << endl;
                }
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
        }

        // Follow the path backward.
        previousEdges.clear();
        if(!isCircular) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                const MarkerGraph::VertexId v0 = edge.source;
                const auto outEdges = markerGraph.edgesBySource[v0];
                if(outEdges.size() != 1) {
                    break;
                }
                const auto inEdges = markerGraph.edgesByTarget[v0];
                if(inEdges.size() != 1) {
                    break;
                }
                edgeId = inEdges[0];
                if(debug) {
                    debugOut << "Backward " << edgeId << " " <<
                        edge.source << "->" << edge.target << endl;
                }
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        if(debug) {
            for(const MarkerGraph::EdgeId edgeId: path) {
                const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                debugOut << "Path " << edgeId << " " <<
                    edge.source << "->" << edge.target << endl;
            }
        }

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            if(wasFound[edgeId]) {
                cout << "Assertion failed at " << edgeId << endl;
                if(debug) {
                    debugOut << "Assertion failed at " << edgeId << endl;
                }
                SHASTA_ASSERT(0);
            }
            wasFound[edgeId] = true;
        }

        // See if this path contains secondary edges.
        bool containsSecondaryEdges = false;
        for(const MarkerGraph::EdgeId edgeId: path) {
            const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
            if(edge.isSecondary) {
                containsSecondaryEdges = true;
                break;
            }
        }

        // Store this path as a new edge of the assembly graph.
        addEdge(path, containsSecondaryEdges);
    }



    // Check that all edges of the marker graph were found.
    SHASTA_ASSERT(find(wasFound.begin(), wasFound.end(), false) == wasFound.end());

}




// Get the vertex descriptor for the vertex corresponding to
// a given MarkerGraph::VertexId, creating the vertex if necessary.
AssemblyGraph2::vertex_descriptor AssemblyGraph2::getVertex(MarkerGraph::VertexId vertexId)
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        const vertex_descriptor v = add_vertex(V(vertexId), *this);
        vertexMap.insert(make_pair(vertexId, v));
        return v;
    } else {
        return it->second;
    }
}



// Create a new edges corresponding to the given path.
// Also create the vertices if necessary.
AssemblyGraph2::edge_descriptor AssemblyGraph2::addEdge(
    const MarkerGraphPath& path,
    bool containsSecondaryEdges)
{
    // Get the first and last edge of the path.
    const MarkerGraph::EdgeId edgeId0 = path.front();
    const MarkerGraph::EdgeId edgeId1 = path.back();
    const MarkerGraph::Edge edge0 = markerGraph.edges[edgeId0];
    const MarkerGraph::Edge edge1 = markerGraph.edges[edgeId1];

    // Get the first and last vertex of the path.
    // This creates the vertices if necessary.
    const MarkerGraph::VertexId vertexId0 = edge0.source;
    const MarkerGraph::VertexId vertexId1 = edge1.target;
    const vertex_descriptor v0 = getVertex(vertexId0);
    const vertex_descriptor v1 = getVertex(vertexId1);

    // Create the edge.
    edge_descriptor e;
    bool edgeWasAdded = false;
    tie(e, edgeWasAdded) = add_edge(v0, v1,
        E(nextEdgeId++, path, containsSecondaryEdges), *this);
    SHASTA_ASSERT(edgeWasAdded);

    return e;
}



void AssemblyGraph2::writeCsv(const string& baseName) const
{
    writeVerticesCsv(baseName + "-Vertices.csv");
    writeEdgesCsv(baseName + "-Edges.csv");
    writeEdgeDetailsCsv(baseName + "-EdgeDetails.csv");
}



void AssemblyGraph2::writeVerticesCsv(const string& fileName) const
{
    const AssemblyGraph2& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0\n";
    BGL_FORALL_VERTICES(v, graph, AssemblyGraph2) {
        csv << graph[v].markerGraphVertexId << "\n";
    }
}



void AssemblyGraph2::writeEdgesCsv(const string& fileName) const
{
    const AssemblyGraph2& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0,VertexId1\n";
    BGL_FORALL_EDGES(e, graph, AssemblyGraph2) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        csv << graph[v0].markerGraphVertexId << ",";
        csv << graph[v1].markerGraphVertexId << "\n";
    }

}



void AssemblyGraph2::writeEdgeDetailsCsv(const string& fileName) const
{
    const G& g = *this;

    ofstream csv(fileName);
    csv << "FirstVertexId,LastVertexId,Branch,Position,EdgeId,VertexId0,VertexId1\n";

    // Loop over edges.
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);
        const MarkerGraph::VertexId vertexId0 = g[v0].markerGraphVertexId;
        const MarkerGraph::VertexId vertexId1 = g[v1].markerGraphVertexId;

        // Loop over branches of this edge.
        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];
            const MarkerGraphPath& path = branch.path;
            for(uint64_t j=0; j<path.size(); j++) {
                const MarkerGraph::EdgeId edgeId = path[j];
                const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                csv << vertexId0 << ",";
                csv << vertexId1 << ",";
                csv << branchId << ",";
                csv << j << ",";
                csv << edgeId << ",";
                csv << edge.source << ",";
                csv << edge.target << "\n";
            }
        }
    }

}





// Assemble sequence for every marker graph path of every edge.
void AssemblyGraph2::assemble()
{
    G& g = *this;

    cout << timestamp << "Assembling sequence." << endl;

    // Use assembled sequence from the marker graph to obtain
    // assembled sequence for all edges.
    BGL_FORALL_EDGES(e, g, G) {
        assemble(e);
    }

    cout << timestamp << "Done assembling sequence." << endl;
}



// Assemble sequence for every marker graph path of a given edge.
void AssemblyGraph2::assemble(edge_descriptor e)
{
    G& g = *this;


    E& edge = g[e];
    for(E::Branch& branch: edge.branches) {
        const MarkerGraphPath& path = branch.path;

        AssembledSegment assembledSegment;
        MarkerGraph::EdgeId const * const begin = &path[0];
        MarkerGraph::EdgeId const * const end = begin + path.size();
        const span<const MarkerGraph::EdgeId> pathSpan(begin, end);
        assembleMarkerGraphPath(k, markers, markerGraph, pathSpan, false, assembledSegment);



        // Store the sequence, excluding the first and last k/2 RLE bases.

        // Compute the number of raw bases to skip at the beginning.
        const uint64_t beginSkip = std::accumulate(
            assembledSegment.repeatCounts.begin(),
            assembledSegment.repeatCounts.begin() + k/2, 0);

        // Compute the number of raw bases to skip at the end.
        const uint64_t endSkip = std::accumulate(
            assembledSegment.repeatCounts.end() - k/2,
            assembledSegment.repeatCounts.end(), 0);

        // Copy the raw bases, excluding those corresponding to the first and last
        // k/2 RLE bases.
        branch.rawSequence.resize(assembledSegment.rawSequence.size() - beginSkip - endSkip);
        copy(
            assembledSegment.rawSequence.begin() + beginSkip,
            assembledSegment.rawSequence.end() - endSkip,
            branch.rawSequence.begin());
    }
}



// Finds edges that form bubbles, then combine
// each of them into a single edge with multiple paths.
void AssemblyGraph2::gatherBubbles()
{
    G& g = *this;



    // Look for sets of parallel edges v0->v1.
    BGL_FORALL_VERTICES(v0, g, G) {

        // Map keyed by child vertex v1, and containing for each v1
        // a vector of all the edges v0->v1.
        std::map< vertex_descriptor, vector<edge_descriptor> > successorMap;
        BGL_FORALL_OUTEDGES(v0, e01, g, G) {
            const vertex_descriptor v1 = target(e01, g);
            successorMap[v1].push_back(e01);
        }



        // Process each set with the same v1.
        for(const auto& p: successorMap) {

            // Get the edges v0->v1.
            const vector<edge_descriptor>& edges01= p.second;
            const uint64_t ploidy = edges01.size();

            // If less than two edges, it is not a bubble, so there is
            // nothing to do.
            if(ploidy < 2) {
                continue;
            }
            const vertex_descriptor v1 = p.first;

            // Create the bubble remove these edges.
            createBubble(v0, v1, edges01);
        }
    }



    // Write out a ploidy histogram.
    vector<uint64_t> ploidyHistogram;
    BGL_FORALL_EDGES(e, g, G) {
        const uint64_t ploidy = g[e].ploidy();
        if(ploidy >= ploidyHistogram.size()) {
            ploidyHistogram.resize(ploidy + 1);
        }
        ++ploidyHistogram[ploidy];
    }
    cout << "Ploidy histogram (counting both strands):" << endl;
    for(uint64_t ploidy=1; ploidy<ploidyHistogram.size(); ploidy++) {
        cout << "Ploidy " << ploidy << ": " << ploidyHistogram[ploidy] << " edges." << endl;
    }

}


AssemblyGraph2::edge_descriptor AssemblyGraph2::createBubble(
    vertex_descriptor v0,
    vertex_descriptor v1,
    const vector<edge_descriptor>& edges01)
{
    G& g = *this;

    // Sanity check.
    for(const edge_descriptor e01: edges01) {
        SHASTA_ASSERT(source(e01, g) == v0);
        SHASTA_ASSERT(target(e01, g) == v1);
        SHASTA_ASSERT(g[e01].ploidy() == 1);
    }

    // Create the new edge to replace the old ones.
    edge_descriptor eNew;
    bool edgeWasAdded = false;
    tie(eNew, edgeWasAdded) = add_edge(v0, v1, E(nextEdgeId++), g);
    SHASTA_ASSERT(edgeWasAdded);
    E& edgeNew = g[eNew];

    // Copy the branches to the new edge and remove the old edges.
    for(const edge_descriptor e01: edges01) {
        const E& edge01 = g[e01];
        copy(edge01.branches.begin(), edge01.branches.end(),
           back_inserter(edgeNew.branches));
        boost::remove_edge(e01, g);
    }

    return eNew;
}



// Find bubbles caused by copy number changes in repeats
// with period up to maxPeriod.
void AssemblyGraph2::findCopyNumberBubbles(uint64_t maxPeriod)
{
    G& g = *this;

    uint64_t totalCount = 0;
    uint64_t bubbleCount = 0;
    uint64_t copyNumberCount = 0;
    BGL_FORALL_EDGES(e, g, G) {
        E& edge = g[e];
        ++totalCount;
        if(not edge.isBubble()) {
            continue;
        }
        ++bubbleCount;
        edge.computeCopyNumberDifferencePeriod(maxPeriod);
        if(edge.period == 0) {
            continue;
        }
        ++copyNumberCount;
        /*
        cout << "Bubble " << edge.id << " of ploidy " << edge.ploidy() <<
            " is a copy number bubble with period " << period << endl;
        */
    }
    cout << "Total number of assembly graph edges " << totalCount << endl;
    cout << "Number of bubbles " << bubbleCount << endl;
    cout << "Number of copy number bubbles " << copyNumberCount << endl;
}



// Double-stranded gfa output.
// This writes a gfa and a csv file with the given base name.
void AssemblyGraph2::writeGfaBothStrands(
    const string& baseName,
    bool writeSequence) const
{
    const G& g = *this;

    // Open the gfa and write the header.
    ofstream gfa(baseName + ".gfa");
    gfa << "H\tVN:Z:1.0\n";

    // Open the csv and write the header.
    ofstream csv(baseName + ".csv");
    csv << ",ComponentId,Phase,Color,First marker graph edge,Last marker graph edge,"
        "Secondary,Period,"
        "Minimum edge coverage,Average edge coverage,Number of distinct oriented reads,\n";



    // Each edge of the AssemblyGraph2 generates a gfa Segment
    // for each of its marker graph paths.
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];

        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];

            // Write a Segment to the GFA file.
            gfa << "S\t" << edge.pathId(branchId) << "\t";
            if(writeSequence) {
                copy(branch.gfaSequence.begin(), branch.gfaSequence.end(),
                    ostream_iterator<Base>(gfa));
                gfa << "\tLN:i:" << branch.gfaSequence.size() << "\n";
            } else {
                gfa << "*\tLN:i:" << branch.path.size() << "\n";
            }

            // Also write a line to the csv file.
            const string color = edge.color(branchId);
            csv <<
                edge.pathId(branchId) << ",";
            if(edge.componentId != std::numeric_limits<uint64_t>::max()) {
                csv << edge.componentId;
            }
            csv << ",";
            if(edge.phase != std::numeric_limits<uint64_t>::max()) {
                csv << (branchId == edge.phase ? 0 : 1);
            }
            csv <<
                "," <<
                color << "," <<
                branch.path.front() << "," << branch.path.back() << "," <<
                (branch.containsSecondaryEdges ? "S" : "") << "," <<
                (edge.period ? to_string(edge.period) : string()) << "," <<
                branch.minimumCoverage << "," <<
                branch.averageCoverage << "," <<
                branch.orientedReadIds.size() << "\n";
        }
    }



    // Generate link record.
    // For each vertex, we generate a Link for each pair of
    // incoming/outgoing marker graph paths.
    BGL_FORALL_VERTICES(v, g, G) {

        // Loop over marker graph paths of incoming edges.
        BGL_FORALL_INEDGES(v, e0, g, G) {
            const E& edge0 = g[e0];
            for(uint64_t i0=0; i0<edge0.ploidy(); i0++) {

                // Loop over marker graph paths of outgoing edges.
                BGL_FORALL_OUTEDGES(v, e1, g, G) {
                    const E& edge1 = g[e1];
                    for(uint64_t i1=0; i1<edge1.ploidy(); i1++) {

                        // To make Bandage happy, we write a Cigar string
                        // consisting of 0M rather than an empty Cigar string.
                        gfa << "L\t" <<
                            edge0.pathId(i0) << "\t+\t" <<
                            edge1.pathId(i1) << "\t+\t0M\n";

                    }
                }
            }
        }
    }
}



string AssemblyGraph2Edge::color(uint64_t branchId) const
{
    if(isBubble()) {

        // A bad or unphased bubble is grey - darker on the strongest branch.
        if(isBad or phase == std::numeric_limits<uint64_t>::max()) {
            if(branchId == strongestBranchId) {
                return "#808080";
            } else {
                return "#C0C0C0";
            }
        }

        if(branchId == phase) {
            return "#bf4040";   // HSL(0,50%,50%) Pink
        } else {
            return "#4040bf";   // HSL(240,50%,50%) Blue
        }

    } else {

        // It is not a bubble.
        return "#808080";
    }

}



// For each bubble edge, compute the number of raw sequence bases
// transferred in each direction for gfa output.
void AssemblyGraph2::countTransferredBases()
{
    G& g = *this;

    BGL_FORALL_EDGES(e, g, G) {
        E& edge = g[e];
        edge.backwardTransferCount = 0;
        edge.forwardTransferCount = 0;

        // To transfer any bases, the edge
        // must be a bubble preceded and followed by a single non-bubble.
        // If these conditions are not satisfied, leave the numbers
        // of transfered bases at 0.

        // The edge must be a bubble.
        if(not edge.isBubble()) {
            continue;
        }

        // v0 must have in-degree and out-degree 1.
        const vertex_descriptor v0 = source(e, g);
        if(in_degree(v0, g) != 1) {
            continue;
        }
        if(out_degree(v0, g) != 1) {
            continue;
        }

        // v1 must have in-degree and out-degree 1.
        const vertex_descriptor v1 = target(e, g);
        if(in_degree(v1, g) != 1) {
            continue;
        }
        if(out_degree(v1, g) != 1) {
            continue;
        }

        // The previous edge must not be a bubble.
        in_edge_iterator itPrevious;
        tie(itPrevious, ignore) = in_edges(v0, g);
        const E& previousEdge = g[*itPrevious];
        if(previousEdge.isBubble()) {
            continue;
        }

        // The next edge must not be a bubble.
        out_edge_iterator itNext;
        tie(itNext, ignore) = out_edges(v1, g);
        const E& nextEdge = g[*itNext];
        if(nextEdge.isBubble()) {
            continue;
        }

        // All conditions are satisfied.
        // Set the number of transfered bases equal to the number
        // of common identical prefix/suffix bases for all the
        // branches of this bubble.
        edge.backwardTransferCount = edge.countCommonPrefixBases();
        edge.forwardTransferCount = edge.countCommonSuffixBases();

        // Make sure we don't transfer more than the length of the
        // shortest branch of this edge.
        uint64_t shortestBranchLength = std::numeric_limits<uint64_t>::max();
        for(const E::Branch& branch:edge.branches) {
            shortestBranchLength = min(shortestBranchLength, uint64_t(branch.rawSequence.size()));
        }
        while(true) {
            if(edge.backwardTransferCount + edge.forwardTransferCount <= shortestBranchLength) {
                break;
            }
            --edge.backwardTransferCount;
            if(edge.backwardTransferCount + edge.forwardTransferCount <= shortestBranchLength) {
                break;
            }
            --edge.forwardTransferCount;
        }
    }
}



// Store GFA sequence in each edge.
// It is obtained from raw sequence by transferring identical bases
// of common sequence between all branches of a bubble to the
// preceding or following non-bubble edge.
void AssemblyGraph2::storeGfaSequence()
{
    G& g = *this;

    // Count the number of sequence bases transferred forward/backward
    // from each bubble edge to adjacent non-bubble edges.
    countTransferredBases();



    BGL_FORALL_EDGES(e, g, G) {
        E& edge = g[e];
        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);

        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            E::Branch& branch = edge.branches[branchId];

            branch.gfaSequence.clear();

            // Add the sequence transferred forward by the preceding bubble, if appropriate.
            if(not edge.isBubble()) {
                if(in_degree(v0, g)==1 and out_degree(v0, g)==1) {
                    in_edge_iterator it;
                    tie(it, ignore) = in_edges(v0, g);
                    const E& previousEdge = g[*it];
                    if(previousEdge.isBubble()) {
                        const vector<Base>& s = previousEdge.branches.front().rawSequence;
                        copy(s.end() - previousEdge.forwardTransferCount, s.end(),
                            back_inserter(branch.gfaSequence));
                    }
                }
            }

            // Add the sequence of this branch, excluding sequence
            // that was transferred backward or forward.
            copy(
                branch.rawSequence.begin() + edge.backwardTransferCount,
                branch.rawSequence.end() - edge.forwardTransferCount,
                back_inserter(branch.gfaSequence));

            // Add the sequence transferred backward by the following bubble, if appropriate.
            if(not edge.isBubble()) {
                if(in_degree(v1, g)==1 and out_degree(v1, g)==1) {
                    out_edge_iterator it;
                    tie(it, ignore) = out_edges(v1, g);
                    const E& nextEdge = g[*it];
                    if(nextEdge.isBubble()) {
                        const vector<Base>& s = nextEdge.branches.front().rawSequence;
                        copy(s.begin(), s.begin() + nextEdge.backwardTransferCount,
                            back_inserter(branch.gfaSequence));
                    }
                }
            }
        }
    }
}



// Return the number of raw bases of sequence identical between
// all branches at the beginning.
uint64_t AssemblyGraph2Edge::countCommonPrefixBases() const
{
    SHASTA_ASSERT(isBubble());

    const vector<Base>& firstRawSequence = branches.front().rawSequence;

    for(uint64_t position=0; position<firstRawSequence.size(); position++) {
        const Base base = firstRawSequence[position];

        for(uint64_t branchId=1; branchId<branches.size(); branchId++) {
            const vector<Base>& rawSequence = branches[branchId].rawSequence;
            if(position == rawSequence.size()) {
                return position;
            }
            if(rawSequence[position] != base) {
                return position;
            }
        }
    }

    return firstRawSequence.size();
}



// Return the number of raw bases of sequence identical between
// all branches at the end.
uint64_t AssemblyGraph2Edge::countCommonSuffixBases() const
{
    SHASTA_ASSERT(isBubble());

    const vector<Base>& firstRawSequence = branches.front().rawSequence;

    for(uint64_t position=0; position<firstRawSequence.size(); position++) {
        const Base base = firstRawSequence[firstRawSequence.size() - 1 - position];

        for(uint64_t branchId=1; branchId<branches.size(); branchId++) {
            const vector<Base>& rawSequence = branches[branchId].rawSequence;
            if(position == rawSequence.size()) {
                return position;
            }
            if(rawSequence[rawSequence.size() - 1 - position] != base) {
                return position;
            }
        }
    }

    return firstRawSequence.size();
}




// Figure out if this is a bubble is caused by copy number
// differences in repeats of period up to maxPeriod.
// If this is the case, returns the shortest period for which this is true.
// Otherwise, returns 0.
void AssemblyGraph2Edge::computeCopyNumberDifferencePeriod(uint64_t maxPeriod)
{
    if(not isBubble()) {
        period = 0;
    }

    // Check all pairs of branches.
    vector<uint64_t> periods;
    for(uint64_t i=0; i<branches.size()-1; i++) {
        const vector<Base>& iSequence = branches[i].rawSequence;
        for(uint64_t j=i+1; j<branches.size(); j++) {
            const vector<Base>& jSequence = branches[j].rawSequence;
            const uint64_t pairPeriod = shasta::isCopyNumberDifference(iSequence, jSequence, maxPeriod);
            if(pairPeriod == 0) {
                period = 0;
                return;
            }
            periods.push_back(pairPeriod);
        }
    }
    deduplicate(periods);


    if(periods.size() == 1) {
        period = periods.front();
    } else {
        period = 0;
    }
}




// Fill in orientedReads and average/minimum coverage.
void AssemblyGraph2Edge::Branch::storeReadInformation(const MarkerGraph& markerGraph)
{
    minimumCoverage = std::numeric_limits<uint64_t>::max();
    averageCoverage = 0;
    orientedReadIds.clear();

    // Loop over the marker graph path of this branch.
    for(MarkerGraph::EdgeId edgeId: path) {

        // Loop over the marker intervals of this edge.
        const span<const MarkerInterval> markerIntervals =
            markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            orientedReadIds.push_back(markerInterval.orientedReadId);
        }

        // Update coverage.
        minimumCoverage = min(minimumCoverage, uint64_t(markerIntervals.size()));
        averageCoverage += markerIntervals.size();
    }

    averageCoverage = uint64_t(std::round(double(averageCoverage) / double(path.size())));
    deduplicate(orientedReadIds);
}



// Store read information on all edges.
void AssemblyGraph2::storeReadInformation()
{
    G& g = *this;

    BGL_FORALL_EDGES(e, g, G) {
        g[e].storeReadInformation(markerGraph);
    }
}


// Store read information on all branches.
void AssemblyGraph2Edge::storeReadInformation(const MarkerGraph& markerGraph)
{
    for(Branch& branch: branches) {
        branch.storeReadInformation(markerGraph);
    }
    findStrongestBranch();
}



void AssemblyGraph2Edge::findStrongestBranch()
{
    SHASTA_ASSERT(not branches.empty());
    strongestBranchId = 0;
    uint64_t strongestBranchCoverage = branches.front().averageCoverage;

    for(uint64_t branchId=0; branchId<branches.size(); branchId++) {
        const uint64_t coverage = branches[branchId].averageCoverage;
        if (coverage > strongestBranchCoverage) {
            strongestBranchId = branchId;
            strongestBranchCoverage = coverage;
        }
    }
}



void AssemblyGraph2::createBubbleGraph(uint64_t readCount)
{
    G& g = *this;

    // Each diploid bubble in the AssemblyGraph2 generates a vertex,
    // except for the repeat count bubbles.
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];

        // If this is not a diploid bubble, skip.
        if(edge.ploidy() != 2) {
            continue;
        }

        // If this is a repeat count bubble, skip.
        if(edge.period != 0) {
            continue;
        }

        // If this is a degenerate bubble, skip.
        if(edge.branches[0].rawSequence == edge.branches[1].rawSequence) {
            continue;
        }

        add_vertex(BubbleGraphVertex(e, edge), bubbleGraph);
    }

    cout << timestamp << "Creating the oriented reads table." << endl;
    bubbleGraph.createOrientedReadsTable(readCount);

    cout << timestamp << "Creating bubble graph edges." << endl;
    bubbleGraph.createEdges();

}



AssemblyGraph2::BubbleGraphVertex::BubbleGraphVertex(
    AssemblyGraph2::edge_descriptor e,
    const AssemblyGraph2Edge& edge) :
    e(e),
    id(edge.id)
{
    SHASTA_ASSERT(edge.ploidy() == 2);
    const vector<OrientedReadId>& orientedReadIds0 = edge.branches[0].orientedReadIds;
    const vector<OrientedReadId>& orientedReadIds1 = edge.branches[1].orientedReadIds;



    // Joint loop over the OrientedReadIds of the two sides.
    const auto begin0 = orientedReadIds0.begin();
    const auto begin1 = orientedReadIds1.begin();
    const auto end0 = orientedReadIds0.end();
    const auto end1 = orientedReadIds1.end();
    auto it0 = begin0;
    auto it1 = begin1;

    while(true) {

        // If both at end, break.
        if(it0 == end0 and it1 == end1) {
            break;
        }

        // If 1 at end but 0 is not, add all remaining orientedReadIds0.
        if(it1==end1) {
            for(; it0 != end0; ++it0) {
                orientedReadIds.push_back(make_pair(*it0, 0));
            }
            break;
        }

        // If 0 at end but 1 is not, add all remaining orientedReadIds1.
        if(it0==end0) {
            for(; it1 != end1; ++it1) {
                orientedReadIds.push_back(make_pair(*it1, 1));
            }
            break;
        }

        // Neither of them is at end.
        if(*it0 < *it1) {
            orientedReadIds.push_back(make_pair(*it0, 0));
            ++it0;
        } else if(*it1 < *it0) {
            orientedReadIds.push_back(make_pair(*it1, 1));
            ++it1;
        } else {
            // It is on both sides. don't store it.
            ++it0;
            ++it1;
        }
    }
}



void AssemblyGraph2::BubbleGraph::createOrientedReadsTable(uint64_t readCount)
{
    BubbleGraph& bubbleGraph = *this;

    orientedReadsTable.clear();
    orientedReadsTable.resize(readCount * 2);
    BGL_FORALL_VERTICES(v, bubbleGraph, BubbleGraph) {
        for(const auto& p: bubbleGraph[v].orientedReadIds) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t side = p.second;
            orientedReadsTable[orientedReadId.getValue()].push_back(make_pair(v, side));
        }
    }
}



void AssemblyGraph2::BubbleGraph::createEdges()
{
    BubbleGraph& bubbleGraph = *this;

    BGL_FORALL_VERTICES(vA, bubbleGraph, BubbleGraph) {
        const BubbleGraphVertex& vertexA = bubbleGraph[vA];
        const uint64_t idA = vertexA.id;

        for(const auto& p: vertexA.orientedReadIds) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t sideA = p.second;

            const auto& v = orientedReadsTable[orientedReadId.getValue()];

            for(const auto& p: v) {
                const vertex_descriptor vB = p.first;
                const BubbleGraphVertex& vertexB = bubbleGraph[vB];
                const uint64_t idB = vertexB.id;

                // Don't add it twice.
                if(idB <= idA) {
                    continue;
                }

                const uint64_t sideB = p.second;

                // Locate the edge between these two bubbles,
                // creating it if necessary.
                BubbleGraph::edge_descriptor e;
                bool edgeWasFound = false;
                tie(e, edgeWasFound) = edge(vA, vB, bubbleGraph);
                if(not edgeWasFound) {
                    tie(e, edgeWasFound) = add_edge(vA, vB, bubbleGraph);
                }
                SHASTA_ASSERT(edgeWasFound);

                // Update the matrix.
                BubbleGraphEdge& edge = bubbleGraph[e];
                ++edge.matrix[sideA][sideB];

            }
        }
    }

}



// Remove edges with too few reads.
void AssemblyGraph2::BubbleGraph::removeWeakEdges(uint64_t minReadCount)
{
    BubbleGraph& bubbleGraph = *this;

    vector<BubbleGraph::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        if(bubbleGraph[e].totalCount() < minReadCount) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const BubbleGraph::edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, bubbleGraph);
    }


}



void AssemblyGraph2::BubbleGraph::computeConnectedComponents()
{
    using boost::connected_components;
    using boost::get;
    using boost::color_map;

    using G = BubbleGraph;
    G& g = *this;

    std::map<vertex_descriptor, uint64_t> colorMap;

    connected_components(g,
        get(&BubbleGraphVertex::componentId, g),
        color_map(boost::make_assoc_property_map(colorMap)));

    // Gather the vertices in each connected component.
    std::map<uint64_t, vector<vertex_descriptor> > componentMap;
    BGL_FORALL_VERTICES(v, g, G) {
        componentMap[g[v].componentId].push_back(v);
    }
    cout << "The BubbleGraph has " << componentMap.size() <<
        " connected components of sizes:";
    connectedComponents.clear();
    for(const auto& p: componentMap) {
        cout << " " << p.second.size();
        connectedComponents.push_back(p.second);
    }
    cout << endl;
}




void AssemblyGraph2::BubbleGraph::writeHtml(const string& fileName)
{
    using G = BubbleGraph;
    G& g = *this;

    /*
    // Create a filtered BubbleGraph, containing only the edges
    // with relativePhase() >= minRelativePhase.
    const double minRelativePhase = 0.;
    using FilteredGraph = boost::filtered_graph<G, BubbleGraphEdgePredicate1>;
    FilteredGraph filteredGraph(g, BubbleGraphEdgePredicate1(g, minRelativePhase));

    // Compute the layout of the filtered graph.
    std::map<FilteredGraph::vertex_descriptor, array<double, 2> > positionMap;
    SHASTA_ASSERT(computeLayout(filteredGraph, "sfdp", 600., positionMap) == ComputeLayoutReturnCode::Success);
    BGL_FORALL_VERTICES(v, filteredGraph, FilteredGraph) {
        filteredGraph[v].position = positionMap[v];
    }
    */

    // Compute the layout.
    std::map<BubbleGraph::vertex_descriptor, array<double, 2> > positionMap;
    SHASTA_ASSERT(computeLayout(g, "sfdp", 600., positionMap) == ComputeLayoutReturnCode::Success);
    BGL_FORALL_VERTICES(v, g, G) {
        g[v].position = positionMap[v];
    }

    // Graphics scaling.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES_T(v, g, G) {
        const auto& position = g[v].position;
        xMin = min(xMin, position[0]);
        xMax = max(xMax, position[0]);
        yMin = min(yMin, position[1]);
        yMax = max(yMax, position[1]);
    }
    const double xyRange = max(xMax-xMin, yMax-yMin);
    const int svgSize = 10000;
    const double vertexRadiusPixels = 3.;
    const double vertexRadius = vertexRadiusPixels * xyRange / double(svgSize);
    const double edgeThicknessPixels = 1.;
    const double edgeThickness = edgeThicknessPixels * xyRange / double(svgSize);

    // Vertex attributes. Color by discordant ratio.
    std::map<G::vertex_descriptor, WriteGraph::VertexAttributes> vertexAttributes;
    BGL_FORALL_VERTICES(v, g, G) {
        const double d = discordantRatio(v);
        const double dRed = 0.3;
        const double hue = max(0., 120. * (1. - d / dRed)); // d=0: green, d=dRed:red
        auto& attributes = vertexAttributes[v];
        attributes.color = "hsl(" + to_string(int(hue)) + ",50%,50%)";
        attributes.radius = vertexRadius;
        attributes.tooltip = to_string(g[v].id);
    }

    // Edge attributes. Color by ambiguity.
    std::map<G::edge_descriptor, WriteGraph::EdgeAttributes> edgeAttributes;
    BGL_FORALL_EDGES(e, g, G) {
        const BubbleGraphEdge& edge = g[e];
        const double ambiguity = edge.ambiguity();
        const double hue = (1. - ambiguity) * 120.; /// Goes from 0 (red) to 120 (green).
        auto& attributes = edgeAttributes[e];
        attributes.color = "hsl(" + to_string(int(hue)) + ",50%,50%)";
        attributes.thickness = edgeThickness;
    }



    // Write it out as svg in html.
    ofstream out(fileName);
    out << "<html><body>";
    WriteGraph::writeSvg(
        g,
        "BubbleGraph",
        svgSize, svgSize,
        vertexAttributes,
        edgeAttributes,
        out);
    out << "</body></html>";
}



void AssemblyGraph2::BubbleGraph::writeGraphviz(const string& fileName) const
{
    const BubbleGraph&  bubbleGraph = *this;

    ofstream out(fileName);
    out <<
        "graph BubbleGraph {\n"
        "node [shape=point];\n";


    // Vertices, colored by discordant ratio.
    // Green if discordant ratio is 0.
    // Red if redDiscordantRatio or more.
    const double redDiscordantRatio = 0.3;
    BGL_FORALL_VERTICES(v, bubbleGraph, BubbleGraph) {
        const BubbleGraphVertex& vertex = bubbleGraph[v];
        const double d = bubbleGraph.discordantRatio(v);
        const double hue = max(0., (1. - d / redDiscordantRatio) / 3.);

        out << vertex.id << " [color=\"" << hue << " 1 1\"];\n";
    }


    // Edges, colored by ambiguity and semi-transparent.
    // Green if ambiguity is 0.
    // Red if ambiguity is redAmbiguity or more.
    const double redAmbiguity = 0.5;
    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        const BubbleGraphEdge& edge = bubbleGraph[e];
        const vertex_descriptor v0 = source(e, bubbleGraph);
        const vertex_descriptor v1 = target(e, bubbleGraph);
        const double ambiguity = edge.ambiguity();
        const double hue = max(0., (1. - ambiguity / redAmbiguity) / 3.);

        out << bubbleGraph[v0].id << "--" << bubbleGraph[v1].id <<
            " [color=\"" << hue << " 1 1 0.5\"];\n";
    }

    out << "}\n";

}



// Post-phasing.
void AssemblyGraph2::BubbleGraph::writeHtml(
    const string& fileName,
    const std::set<edge_descriptor>& treeEdges,
    bool onlyShowUnhappyEdges)
{
    BubbleGraph& g = *this;

    // Create a filtered BubbleGraph, containing only the tree edges
    // and the edges with relativePhase() >= minRelativePhase.
    const double minRelativePhase = 0.;
    using FilteredGraph = boost::filtered_graph<BubbleGraph, BubbleGraphEdgePredicate2>;
    FilteredGraph filteredGraph(g, BubbleGraphEdgePredicate2(g, minRelativePhase, treeEdges));

    // Compute the layout of the filtered graph.
    // This should mostly separate the two phases.
    std::map<FilteredGraph::vertex_descriptor, array<double, 2> > positionMap;
    SHASTA_ASSERT(computeLayout(filteredGraph, "sfdp", 600., positionMap) == ComputeLayoutReturnCode::Success);
    BGL_FORALL_VERTICES(v, filteredGraph, FilteredGraph) {
        filteredGraph[v].position = positionMap[v];
    }



    // Graphics scaling.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES_T(v, g, G) {
        const auto& position = g[v].position;
        xMin = min(xMin, position[0]);
        xMax = max(xMax, position[0]);
        yMin = min(yMin, position[1]);
        yMax = max(yMax, position[1]);
    }
    const double xyRange = max(xMax-xMin, yMax-yMin);
    const int svgSize = 1600;
    const double vertexRadiusPixels = 3.;
    const double vertexRadius = vertexRadiusPixels * xyRange / double(svgSize);
    const double edgeThicknessPixels = 0.3;
    const double edgeThickness = edgeThicknessPixels * xyRange / double(svgSize);

    // Vertex attributes. Color by phase.
    std::map<G::vertex_descriptor, WriteGraph::VertexAttributes> vertexAttributes;
    BGL_FORALL_VERTICES(v, g, BubbleGraph) {
        auto& attributes = vertexAttributes[v];
        if(g[v].phase == 0) {
            attributes.color = "hsl(240,50%,50%)";
        } else {
            attributes.color = "hsl(300,50%,50%)";
        }
        attributes.radius = vertexRadius;
        attributes.tooltip = to_string(g[v].id);
    }

    // Edge attributes.
    std::map<BubbleGraph::edge_descriptor, WriteGraph::EdgeAttributes> edgeAttributes;
    BGL_FORALL_EDGES(e, g, BubbleGraph) {
        const double relativePhase = g[e].relativePhase();
        const double hue = (1. + relativePhase) * 60.; /// Goes from 0 (red) to 120 (green).
        auto& attributes = edgeAttributes[e];
        attributes.color = "hsla(" + to_string(int(hue)) + ",50%,50%,50%)";
        attributes.thickness = edgeThickness;
        if(onlyShowUnhappyEdges) {
            if(g.edgeIsHappy(e)) {
                attributes.color = "hsla(" + to_string(int(hue)) + ",50%,50%,0%)";
            } else {
                attributes.thickness *= 20.;
            }
        }
    }



    // Write it out as svg in html.
    ofstream out(fileName);
    out << "<html><body>";
    WriteGraph::writeSvg(
        g,
        "BubbleGraph",
        svgSize, svgSize,
        vertexAttributes,
        edgeAttributes,
        out);
    out << "</body></html>";
}



// Return true if the give edge has relative phase consistent
// with the phases assigned to its two vertices.
bool AssemblyGraph2::BubbleGraph::edgeIsHappy(BubbleGraph::edge_descriptor e) const
{
    const BubbleGraph& g = *this;

    const vertex_descriptor v0 = source(e, g);
    const vertex_descriptor v1 = target(e, g);
    const uint64_t phase0 = g[v0].phase;
    const uint64_t phase1 = g[v1].phase;

    const double relativePhase = g[e].relativePhase();

    if(phase0 == phase1) {
        return relativePhase > 0.;
    } else {
        return relativePhase < 0.;
    }
}



double AssemblyGraph2::BubbleGraph::discordantRatio(vertex_descriptor v) const
{
    const BubbleGraph& bubbleGraph = *this;

    uint64_t concordantSum = 0;
    uint64_t discordantSum = 0;
    BGL_FORALL_OUTEDGES(v, e, bubbleGraph, BubbleGraph) {
        const BubbleGraphEdge& edge = bubbleGraph[e];
        concordantSum += edge.concordantCount();
        discordantSum += edge.discordantCount();
    }

    return double(discordantSum) / double(concordantSum + discordantSum);
}



void AssemblyGraph2::BubbleGraph::removeWeakVertices(
    double discordantRatioThreshold,
    vector<AssemblyGraph2::edge_descriptor>& badBubbles)
{
    BubbleGraph& bubbleGraph = *this;

    vector<BubbleGraph::vertex_descriptor> verticesToBeRemoved;
    badBubbles.clear();
    BGL_FORALL_VERTICES(v, bubbleGraph, BubbleGraph) {
        if(discordantRatio(v) > discordantRatioThreshold) {
            verticesToBeRemoved.push_back(v);
            badBubbles.push_back(bubbleGraph[v].e);
        }
    }

    for(const BubbleGraph::vertex_descriptor v: verticesToBeRemoved) {
        clear_vertex(v, bubbleGraph);
        remove_vertex(v, bubbleGraph);
    }
}



void AssemblyGraph2::BubbleGraph::writeEdgesCsv(const string& fileName) const
{
    const BubbleGraph& bubbleGraph = *this;

    ofstream csv(fileName);
    csv << "BubbleIdA,BubbleIdB,m00,m11,m01,m10\n";
    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        const BubbleGraphEdge& edge = bubbleGraph[e];
        uint64_t idA = bubbleGraph[source(e, bubbleGraph)].id;
        uint64_t idB = bubbleGraph[target(e, bubbleGraph)].id;
        if(idB < idA) {
            std::swap(idA, idB);
        }
        csv << idA << ",";
        csv << idB << ",";
        csv << edge.matrix[0][0] << ",";
        csv << edge.matrix[1][1] << ",";
        csv << edge.matrix[0][1] << ",";
        csv << edge.matrix[1][0] << "\n";
    }

}



// Create a new BubbleGraph from a given connected component.
void AssemblyGraph2::BubbleGraph::extractComponent(
    uint64_t componentId,
    BubbleGraph& componentGraph) const
{
    const BubbleGraph& bubbleGraph = *this;

    SHASTA_ASSERT(componentId < connectedComponents.size());
    const vector<vertex_descriptor>& componentVertices = connectedComponents[componentId];

    SHASTA_ASSERT(num_vertices(componentGraph) == 0);
    SHASTA_ASSERT(num_edges(componentGraph) == 0);

    // Map with:
    // Key = vertex_descriptor in the full graph.
    // Value = vertex descriptor in the component graph.
    std::map<vertex_descriptor, vertex_descriptor> vertexMap;

    // Create the vertices.
    for(const vertex_descriptor v: componentVertices) {
        const vertex_descriptor vc = add_vertex(bubbleGraph[v], componentGraph);
        vertexMap.insert(make_pair(v, vc));
    }

    // Create the edges.
    for(const vertex_descriptor v0: componentVertices) {
        const auto it0 = vertexMap.find(v0);
        SHASTA_ASSERT(it0 != vertexMap.end());
        const vertex_descriptor vc0 = it0->second;
        BGL_FORALL_OUTEDGES(v0, e01, bubbleGraph, BubbleGraph) {
            const vertex_descriptor v1 = target(e01, bubbleGraph);
            if(bubbleGraph[v0].id < bubbleGraph[v1].id) {
                const auto it1 = vertexMap.find(v1);
                SHASTA_ASSERT(it1 != vertexMap.end());
                const vertex_descriptor vc1 = it1->second;
                add_edge(vc0, vc1, bubbleGraph[e01], componentGraph);
            }
        }
    }
}




// Use each connected component of the bubble graph to phase the bubbles.
void AssemblyGraph2::phase()
{
    const bool debug = false;

    for(uint64_t componentId=0;
        componentId<bubbleGraph.connectedComponents.size(); componentId++) {
        BubbleGraph componentGraph;
        bubbleGraph.extractComponent(componentId, componentGraph);

        cout << "Processing connected component " << componentId <<
            " with " << num_vertices(componentGraph) <<
            " vertices and " << num_edges(componentGraph) <<
            " edges." << endl;

        if(debug) {
            componentGraph.writeGraphviz("Component-" + to_string(componentId) + ".dot");
        }

        // Compute an index map, needed below, which maps vertices to integers.
        std::map<BubbleGraph::vertex_descriptor, uint64_t> indexMap;
        uint64_t vertexIndex = 0;
        BGL_FORALL_VERTICES(v, componentGraph, BubbleGraph) {
            indexMap.insert(make_pair(v, vertexIndex++));
        }

        // Compute a minimal spanning tree that minimizes
        // the sum of  edge weights defined as
        // difference discordantCount() - concordantCount()
        std::map<BubbleGraph::edge_descriptor, int64_t> weightMap;
        BGL_FORALL_EDGES(e, componentGraph, BubbleGraph) {
            const BubbleGraphEdge& edge = componentGraph[e];
            const int64_t weight = int64_t(edge.discordantCount()) - int64_t(edge.concordantCount());
            weightMap.insert(make_pair(e, weight));
        }
        std::set<BubbleGraph::edge_descriptor> treeEdges;
        boost::kruskal_minimum_spanning_tree(componentGraph, std::inserter(treeEdges, treeEdges.begin()),
            weight_map(boost::make_assoc_property_map(weightMap)).
            vertex_index_map(boost::make_assoc_property_map(indexMap)));
        SHASTA_ASSERT(treeEdges.size() == indexMap.size() - 1);



        // Write out the tree edges to csv.
        if(debug) {
            ofstream csv("Component-" + to_string(componentId) + "-TreeEdges.csv");
            csv << "BubbleId0,BubbleId1,Diagonal,OffDiagonal,Concordant,Discordant,Weight\n";
            for(const BubbleGraph::edge_descriptor e: treeEdges) {
                const BubbleGraphEdge& edge = componentGraph[e];
                const BubbleGraph::vertex_descriptor v0 = source(e, bubbleGraph);
                const BubbleGraph::vertex_descriptor v1 = target(e, bubbleGraph);
                csv << bubbleGraph[v0].id << ",";
                csv << bubbleGraph[v1].id << ",";
                csv << edge.diagonalCount() << ",";
                csv << edge.offDiagonalCount() << ",";
                csv << edge.concordantCount() << ",";
                csv << edge.discordantCount() << ",";
                csv << edge.concordantCount() - edge.discordantCount() << "\n";
            }
        }



        // To phase, do a BFS on the spanning tree of the componentGraph.
        // Assign phases consistent with the spanning tree edges.
        std::queue<BubbleGraph::vertex_descriptor> q;
        BubbleGraph::vertex_iterator it;
        tie(it, ignore) = vertices(componentGraph);
        BubbleGraph::vertex_descriptor vStart = *it;
        q.push(vStart);
        componentGraph[vStart].phase = 0;
        while(not q.empty()) {

            // Dequeue a vertex.
            const BubbleGraph::vertex_descriptor v0 = q.front();
            // cout << "Dequeued " << (*this)[componentGraph[v0].e].id << " " << componentGraph[v0].phase << endl;
            q.pop();
            const uint64_t phase0 = componentGraph[v0].phase;
            SHASTA_ASSERT(phase0 != BubbleGraphVertex::invalidPhase);

            // Loop over tree edges.
            BGL_FORALL_OUTEDGES(v0, e01, componentGraph, BubbleGraph) {
                if(treeEdges.find(e01) == treeEdges.end()) {
                    continue;
                }
                const double relativePhase = componentGraph[e01].relativePhase();
                const BubbleGraph::vertex_descriptor v1 = target(e01, componentGraph);
                // cout << "Found " << (*this)[componentGraph[v1].e].id << " " << componentGraph[v1].phase << endl;
                // cout << "Relative phase " << relativePhase << endl;
                uint64_t& phase1 = componentGraph[v1].phase;
                if(phase1 == BubbleGraphVertex::invalidPhase) {
                    if(relativePhase > 0.) {
                        phase1 = phase0;
                    } else {
                        phase1 = 1 - phase0;
                    }
                    q.push(v1);
                    // cout << "Enqueued " << (*this)[componentGraph[v1].e].id << " " << componentGraph[v1].phase << endl;
                } else {
                    // We already encountered this vertex. Just check
                    // that its phase is consistent with the edge.
                    const uint64_t& phase1 = componentGraph[v1].phase;
                    if(relativePhase > 0.) {
                        SHASTA_ASSERT(phase1 == phase0);
                    } else {
                        SHASTA_ASSERT(phase1 == 1 - phase0);
                    }
                }
            }
        }

        uint64_t unhappyCount = 0;
        uint64_t totalCount = 0;
        BGL_FORALL_EDGES(e, componentGraph, BubbleGraph) {
            ++totalCount;
            if(not componentGraph.edgeIsHappy(e)) {
                ++unhappyCount;
            }
        }
        cout << "Found " << unhappyCount << " edges inconsistent with computed bubble phasing "
            " out of " << totalCount << " edges in this connected component." << endl;

        if(debug) {
            componentGraph.writeHtml("Component-" + to_string(componentId) + ".html", treeEdges, false);
            componentGraph.writeHtml("Component-" + to_string(componentId) + "-Unhappy.html", treeEdges, true);
        }


        // Copy the phasing to the global bubble graph.
        uint64_t i = 0;
        BGL_FORALL_VERTICES(v, componentGraph, BubbleGraph) {
            BubbleGraph::vertex_descriptor vGlobal = bubbleGraph.connectedComponents[componentId][i++];
            bubbleGraph[vGlobal].phase = bubbleGraph[v].phase;
            SHASTA_ASSERT(bubbleGraph[vGlobal].componentId == componentId);
        }
    }


    // Check that all vertices of the global bubble graph are phased
    // and copy the phasing information to the AssemvblyGraph2Edges.
    BGL_FORALL_VERTICES(v, bubbleGraph, BubbleGraph) {
        const BubbleGraphVertex& bubbleGraphVertex = bubbleGraph[v];
        const uint64_t componentId = bubbleGraphVertex.componentId;
        SHASTA_ASSERT(componentId != std::numeric_limits<uint64_t>::max());
        const uint64_t phase = bubbleGraphVertex.phase;
        SHASTA_ASSERT(phase != BubbleGraphVertex::invalidPhase);

        AssemblyGraph2Edge& assemblyGraph2Edge = (*this)[bubbleGraphVertex.e];
        assemblyGraph2Edge.componentId = componentId;
        assemblyGraph2Edge.phase = phase;
    }
}



void AssemblyGraph2::removeSecondaryBubbles()
{
    G& g = *this;

    // Loop over edges.
    BGL_FORALL_EDGES(e, g, G) {

        // If not a bubble, do nothing.
        E& edge = g[e];
        if(not edge.isBubble()) {
            continue;
        }

        // See how many secondary branches we have.
        uint64_t secondaryCount = 0;
        for(const E::Branch& branch: edge.branches) {
            if(branch.containsSecondaryEdges) {
                ++secondaryCount;
            }
        }

        // If there are no secondary branches, do nothing.
        if(secondaryCount == 0) {
            continue;
        }
        const uint64_t primaryCount = edge.ploidy() - secondaryCount;


        // Remove secondary branches, keeping at most one.
        if(primaryCount > 0) {

            // There is at least one primary branch, so we can
            // remove all the secondary ones.
            edge.removeAllSecondaryBranches();
            edge.findStrongestBranch();
        } else {

            // There are no primary branches.
            // Remove all secondary branches except the strongest.
            edge.removeAllBranchesExceptStrongest();
            edge.findStrongestBranch();
        }
    }

}



void AssemblyGraph2Edge::removeAllSecondaryBranches()
{
    vector<Branch> newBranches;
    for(const Branch& branch: branches) {
        if(not branch.containsSecondaryEdges) {
            newBranches.push_back(branch);
        }
    }
    branches.swap(newBranches);
}



void AssemblyGraph2Edge::removeAllBranchesExceptStrongest()
{
    vector<Branch> newBranches(1, branches[strongestBranchId]);
    branches.swap(newBranches);
}



// Merge consecutive non-bubbles, when possible.
void AssemblyGraph2::merge()
{
    // Find linear chains of non-bubbles.
    vector< vector<edge_descriptor> > chains;
    findNonBubbleLinearChains(chains);

    /*
    cout << "Found the following linear chains that will be merged:" << endl;
    for(const auto& chain: chains) {
        cout << "Chain";
        for(const edge_descriptor e: chain) {
            cout << " " << (*this)[e].id;
        }
        cout << endl;
    }
    */

    // Merge each chain.
    for(const vector<edge_descriptor>& chain: chains) {
        merge(chain);
    }
}



AssemblyGraph2::edge_descriptor AssemblyGraph2::merge(const vector<edge_descriptor>& chain)
{
    G& g = *this;

    /*
    cout << "Merging linear chain ";
    for(const edge_descriptor e: chain) {
        cout << " " << g[e].id;
    }
    cout << endl;
    */


    // Sanity checks.

    // Must be more than one edge.
    SHASTA_ASSERT(chain.size() > 1);

    // There must be no bubbles.
    for(const edge_descriptor e: chain) {
        SHASTA_ASSERT(not g[e].isBubble());
    }

    // It must be a path.
    for(uint64_t i=1; i<chain.size(); i++) {
        SHASTA_ASSERT(target(chain[i-1], g) == source(chain[i], g));
    }

    // Check the degrees of internal vertices.
    // It must be a linear chain.
    for(uint64_t i=1; i<chain.size(); i++) {
        const vertex_descriptor v = source(chain[i], g);
        SHASTA_ASSERT(in_degree(v, g) == 1);
        SHASTA_ASSERT(out_degree(v, g) == 1);
    }

    // Create the merged marker graph path.
    MarkerGraphPath newPath;
    bool containsSecondaryEdges = false;
    for(const edge_descriptor e: chain) {
        const AssemblyGraph2Edge::Branch& branch = g[e].branches.front();
        const MarkerGraphPath& path = g[e].branches.front().path;
        newPath.insert(newPath.end(), path.begin(), path.end());
        containsSecondaryEdges = containsSecondaryEdges or branch.containsSecondaryEdges;
    }

    // Add the new edge.
    const edge_descriptor eNew = addEdge(newPath, containsSecondaryEdges);
    g[eNew].storeReadInformation(markerGraph);

    // Gather the vertices in between, which will be removed.
    vector<vertex_descriptor> verticesToBeRemoved;
    for(uint64_t i=1; i<chain.size(); i++) {
        verticesToBeRemoved.push_back(source(chain[i], g));
    }

    // Remove the old edges.
    for(const edge_descriptor e: chain) {
        boost::remove_edge(e, g);
    }

    // Remove the vertices in between.
    for(const vertex_descriptor v: verticesToBeRemoved) {
        SHASTA_ASSERT(in_degree(v, g) == 0);
        SHASTA_ASSERT(out_degree(v, g) == 0);
        boost::remove_vertex(v, g);
    }

    return eNew;
}



// Find linear chains of adjacent non-bubbles.
// Used by merge.
void AssemblyGraph2::findNonBubbleLinearChains(
    vector< vector<edge_descriptor> >& chains) const
{
    const G& g = *this;

    // The edges we have already encountered.
    std::set<edge_descriptor> edgesFound;

    // Start with no chains.
    chains.clear();

    // Work vectors used in the main loop below.
    vector<edge_descriptor> forwardPath;
    vector<edge_descriptor> backwardPath;

    // Consider all possible start edges for the chain.
    BGL_FORALL_EDGES_T(eStart, g, G) {

        // If this is a bubble, skip it.
        if(g[eStart].isBubble()) {
            continue;
        }

        // If we already assigned this edge to a chain, skip it.
        if(edgesFound.find(eStart) != edgesFound.end()) {
            continue;
        }

        edgesFound.insert(eStart);

        // Extend forward.
        forwardPath.clear();
        bool isCircular = false;
        edge_descriptor e = eStart;
        while(true) {
            const vertex_descriptor v = target(e, g);
            if(in_degree(v, g) != 1) {
                break;
            }
            if(out_degree(v, g) != 1) {
                break;
            }
            BGL_FORALL_OUTEDGES_T(v, eNext, g, G) {
                e = eNext;
                break;
            }
            if(g[e].isBubble()) {
                break;
            }
            if(e == eStart) {
                isCircular = true;
                break;
            }
            forwardPath.push_back(e);
            SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
            edgesFound.insert(e);
        }


        // Extend backward.
        backwardPath.clear();
        if(not isCircular) {
            edge_descriptor e = eStart;
            while(true) {
                const vertex_descriptor v = source(e, g);
                if(in_degree(v, g) != 1) {
                    break;
                }
                if(out_degree(v, g) != 1) {
                    break;
                }
                BGL_FORALL_INEDGES_T(v, ePrevious, g, G) {
                    e = ePrevious;
                    break;
                }
                if(g[e].isBubble()) {
                    break;
                }
                if(e == eStart) {
                    isCircular = true;
                    break;
                }
                backwardPath.push_back(e);
                SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
                edgesFound.insert(e);
            }
        }

        // If the chain is too short, get rid of it.
        if(backwardPath.size() + 1 + forwardPath.size() < 2) {
            continue;
        }

        // Store it.
        chains.resize(chains.size() + 1);
        vector<edge_descriptor>& chain = chains.back();
        copy(backwardPath.rbegin(), backwardPath.rend(), back_inserter(chain));
        chain.push_back(eStart);
        copy(forwardPath.begin(), forwardPath.end(), back_inserter(chain));

    }

    // Check that all non-bubble edges were found.
    BGL_FORALL_EDGES_T(e, g, G) {
        if(not g[e].isBubble()) {
            SHASTA_ASSERT(edgesFound.find(e) != edgesFound.end());
        }
    }

}


