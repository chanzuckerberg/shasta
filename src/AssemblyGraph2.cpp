#include "AssemblyGraph2.hpp"
#include "AssembledSegment.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
using namespace shasta;

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

    // Maximum period to remove bubbles caused by copy number differences.
    const uint64_t maxPeriod = 4;

    // Create the AssemblyGraph2 from the MarkerGraph,
    // gather each set of bubbles in a single edge,
    // and assemble sequence on every edge.
    create();
    gatherBubbles();
    assemble();
    storeGfaSequence();
    checkReverseComplementEdges();

    // Find bubbles caused by copy number changes in repeats
    // with period up to maxPeriod.
    findCopyNumberBubbles(maxPeriod);

    // Store read information on all edges.
    storeReadInformation();

    // Write out what we have.
    const bool writeSequence = true;
    writeGfaBothStrands("Assembly-BothStrands", writeSequence);
    writeGfa("Assembly", writeSequence);

    // Diploid phasing of the bubbles.
    createBubbleGraph(markers.size()/2);


}



// Initial creation of vertices and edges.
void AssemblyGraph2::create()
{
    G& g = *this;

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
                startEdge.source << "->" << startEdge.target << "\n";
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
            if(edgeId == startEdgeId) {
                isCircular = true;
                if(true) {
                    cout << "Found a circular edge." << endl;
                }
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
            if(debug) {
                debugOut << "Forward " << edgeId << " " <<
                    edge.source << "->" << edge.target << "\n";

            }
        }

        // Follow the path backward.
        if(!isCircular) {
            previousEdges.clear();
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
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
                if(debug) {
                    debugOut << "Backward " << edgeId << " " <<
                        edge.source << "->" << edge.target << "\n";
                }
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            wasFound[edgeId] = true;
        }

        // Store this path as a new edge of the assembly graph.
        const edge_descriptor e = addEdge(path);
        if(debug) {
            for(const MarkerGraph::EdgeId edgeId: path) {
                const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                debugOut << "Path " << edgeId << " " <<
                    edge.source << "->" << edge.target << "\n";
            }
        }

        // Also construct the reverse complemented path.
        reverseComplementedPath.clear();
        for(const MarkerGraph::EdgeId edgeId: path) {
            const MarkerGraph::EdgeId edgeIdRc = markerGraph.reverseComplementEdge[edgeId];
            reverseComplementedPath.push_back(edgeIdRc);
        }
        std::reverse(reverseComplementedPath.begin(), reverseComplementedPath.end());



        // Figure out if the reverse complemented chain is the same
        // as the original chain. This can happen in exceptional cases.
        bool isSelfComplementary = false;
        if(!isCircular) {
            isSelfComplementary = (path == reverseComplementedPath);
        } else {

            // For a circular path the test is more complex.
            // We check if the reverse complement of the first edge
            // is in the path.
            isSelfComplementary =
                find(path.begin(), path.end(), reverseComplementedPath.front()) != path.end();
        }


        // Store the reverse complemented path, if different from the original one.
        if(not isSelfComplementary) {
            const edge_descriptor eRc = addEdge(reverseComplementedPath);
            for(const MarkerGraph::EdgeId edgeIdRc: reverseComplementedPath) {
                SHASTA_ASSERT(not wasFound[edgeIdRc]);
                wasFound[edgeIdRc] = true;
            }

            // The two edges we added are reverse complement of each other.
            g[e].reverseComplement = eRc;
            g[eRc].reverseComplement = e;

            if(debug) {
                for(const MarkerGraph::EdgeId edgeId: path) {
                    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                    debugOut << "Reverse complemented path " << edgeId << " " <<
                        edge.source << "->" << edge.target << "\n";
                }
            }
        } else {
            // The edge we added is reverse complement of itself.
            g[e].reverseComplement = e;
            cout << "Found a self-complementary edge." << endl;
        }

    }



    // Check that all edges of the marker graph were found.
    SHASTA_ASSERT(find(wasFound.begin(), wasFound.end(), false) == wasFound.end());

}



void AssemblyGraph2::checkReverseComplementEdges() const
{
    const G& g = *this;

    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        const edge_descriptor eRc = edge.reverseComplement;
        const E& edgeRc = g[eRc];
        SHASTA_ASSERT(edgeRc.reverseComplement == e);

        if(e != eRc) {

            SHASTA_ASSERT(edge.ploidy() == edgeRc.ploidy());
            for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
                const E::Branch& branch = edge.branches[branchId];
                const E::Branch& branchRc = edgeRc.branches[branchId];

                // Check that the paths are the reverse complement of each other.
                // For a circular path the check is more complicated, so don't bother.
                const MarkerGraphPath& path = branch.path;

                const bool isCircular =
                    markerGraph.edges[path.front()].source ==
                    markerGraph.edges[path.back()].target;

                if(not isCircular) {
                    const MarkerGraphPath& pathRc = branchRc.path;
                    const uint64_t n = path.size();
                    SHASTA_ASSERT(pathRc.size() == n);
                    for(uint64_t i=0; i<n; i++) {
                        const MarkerGraph::EdgeId edgeId = path[i];
                        const MarkerGraph::EdgeId edgeIdRc = pathRc[n - 1 - i];
                        SHASTA_ASSERT(markerGraph.reverseComplementEdge[edgeId] == edgeIdRc);
                        SHASTA_ASSERT(markerGraph.reverseComplementEdge[edgeIdRc] == edgeId);
                    }

                    // Also check that the sequences are reverse complement of each other.
                    const vector<Base>& sequence = branch.rawSequence;
                    const vector<Base>& sequenceRc = branchRc.rawSequence;
                    SHASTA_ASSERT(sequence.size() == sequenceRc.size());
                    for(uint64_t position=0; position<sequence.size(); position++) {
                        const Base base = sequence[position];
                        const Base baseRc = sequenceRc[sequence.size() - 1 - position];
                        SHASTA_ASSERT(baseRc == base.complement());
                    }

                    // Also check that the gfa sequences are reverse complement of each other.
                    const vector<Base>& gfaSequence = branch.gfaSequence;
                    const vector<Base>& gfaSequenceRc = branchRc.gfaSequence;
                    SHASTA_ASSERT(gfaSequence.size() == gfaSequenceRc.size());
                    for(uint64_t position=0; position<gfaSequence.size(); position++) {
                        const Base base = gfaSequence[position];
                        const Base baseRc = gfaSequenceRc[gfaSequence.size() - 1 - position];
                        SHASTA_ASSERT(baseRc == base.complement());
                    }
                }
            }
        }
    }
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
AssemblyGraph2::edge_descriptor AssemblyGraph2::addEdge(const MarkerGraphPath& path)
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
    tie(e, edgeWasAdded) = add_edge(v0, v1, E(nextEdgeId++, path), *this);
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
// This guarantees that, for each pair of reverse complemented edges,
// The assembled sequences are the reverse complement of each other.
void AssemblyGraph2::assemble()
{
    G& g = *this;

    cout << timestamp << "Assembling sequence." << endl;

    // Use assembled sequence from the marker graph to obtain
    // assembled sequence for all edges that have an id
    // not greater than the id of their reverse complement.
    BGL_FORALL_EDGES(e, g, G) {
        if(not idIsGreaterThanReverseComplement(e)) {
            assemble(e);
        }
    }



    // For the remaining edges, obtain the sequence
    // from their reverse complement.
    BGL_FORALL_EDGES(e, g, G) {
        if(idIsGreaterThanReverseComplement(e)) {

            // Access this edge and its reverse complement.
            E& edge = g[e];
            const edge_descriptor eRc = edge.reverseComplement;
            const E& edgeRc = g[eRc];

            // Sanmity check: they must have the same ploidy.
            const uint64_t ploidy = edge.ploidy();
            SHASTA_ASSERT(edgeRc.ploidy() == ploidy);

            // For each branch, copy the sequence from the
            // reverse complement edge, while reverse complementing.
            for(uint64_t branchId=0; branchId<ploidy; branchId++) {
                const vector<Base>& sequenceRc = edgeRc.branches[branchId].rawSequence;
                vector<Base>& sequence = edge.branches[branchId].rawSequence;
                sequence.resize(sequenceRc.size());
                copy(sequenceRc.rbegin(), sequenceRc.rend(), sequence.begin());
                for(Base& b: sequence) {
                    b.complementInPlace();
                }
            }
        }
    }

    cout << timestamp << "Done assembling sequence." << endl;
}



// Return true if an edge has id less than its reverse complement.
bool AssemblyGraph2::idIsLessThanReverseComplement(edge_descriptor e) const
{
    const G& g = *this;

    const E& edge = g[e];
    const edge_descriptor eRc = edge.reverseComplement;
    const E& edgeRc = g[eRc];
    return edge.id < edgeRc.id;
}



// Return true if an edge has id greater than its reverse complement.
bool AssemblyGraph2::idIsGreaterThanReverseComplement(edge_descriptor e) const
{
    const G& g = *this;

    const E& edge = g[e];
    const edge_descriptor eRc = edge.reverseComplement;
    const E& edgeRc = g[eRc];
    return edge.id > edgeRc.id;
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

            // Also gather the reverse complement edges, in the same order.
            vector<edge_descriptor> edges01Rc(ploidy);
            for(uint64_t branchId=0; branchId<ploidy; branchId++) {
                const edge_descriptor e01 = edges01[branchId];
                const edge_descriptor e01Rc = g[e01].reverseComplement;
                edges01Rc[branchId] = e01Rc;
            }

            // Get the vertices of the reverse complemented edges.
            const vertex_descriptor v1Rc = source(edges01Rc.front(), g);
            const vertex_descriptor v0Rc = target(edges01Rc.front(), g);
            for(const edge_descriptor e01Rc: edges01Rc) {
                SHASTA_ASSERT(source(e01Rc, g) == v1Rc);
                SHASTA_ASSERT(target(e01Rc, g) == v0Rc);
            }


            // Create the bubble and its reverse complement and remove these edges.
            const edge_descriptor eNew = createBubble(v0, v1, edges01);
            if(v0Rc == v0) {
                // The bubble is the same as its reverse complement.
                cout << "Found a self-complementary bubble." << endl;
                g[eNew].reverseComplement = eNew;
            } else {
                const edge_descriptor eNewRc = createBubble(v1Rc, v0Rc, edges01Rc);
                g[eNew].reverseComplement = eNewRc;
                g[eNewRc].reverseComplement = eNew;

            }
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
    csv << ",Color,First marker graph edge,Last marker graph edge,"
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
            string color = "Grey";
            if(edge.isBubble()) {
                color = edge.colorByPeriod(branchId);
            }
            csv <<
                edge.pathId(branchId) << "," <<
                color << "," <<
                branch.path.front() << "," << branch.path.back() << "," <<
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



// Single-stranded gfa output.
// This writes a gfa and a csv file with the given base name.
void AssemblyGraph2::writeGfa(
    const string& baseName,
    bool writeSequence) const
{
    const G& g = *this;

    // Open the gfa and write the header.
    ofstream gfa(baseName + ".gfa");
    gfa << "H\tVN:Z:1.0\n";

    // Open the csv and write the header.
    ofstream csv(baseName + ".csv");
    csv << ",Color,First marker graph edge,Last marker graph edge,"
        "Minimum edge coverage,Average edge coverage,Number of distinct oriented reads,\n";



    // Each edge of the AssemblyGraph2 generates a gfa Segment
    // for each of its marker graph paths.
    // For single strand output, skip edges with id
    // greater than the id of their reverse complement.
    BGL_FORALL_EDGES(e, g, G) {
        if(idIsGreaterThanReverseComplement(e)) {
            continue;
        }
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
            string color = "Grey";
            if(edge.isBubble()) {
                color = edge.colorByPeriod(branchId);
            }
            csv <<
                edge.pathId(branchId) << "," <<
                color << "," <<
                branch.path.front() << "," << branch.path.back() << "," <<
                branch.minimumCoverage << "," <<
                branch.averageCoverage << "," <<
                branch.orientedReadIds.size() << "\n";
        }
    }


    // Generate link record.
    // For each vertex, we generate a Link for each pair of
    // incoming/outgoing marker graph paths.
    // Skip as necessary to avoid writing links twice.
    BGL_FORALL_VERTICES(v, g, G) {

        // Loop over branches of incoming edges.
        BGL_FORALL_INEDGES(v, e0, g, G) {

            // Reverse complement e0 if necessary.
            edge_descriptor ee0 = e0;
            bool reverse0 = false;
            if(idIsGreaterThanReverseComplement(e0)) {
                ee0 = g[e0].reverseComplement;
                reverse0 = true;
            }
            const E& edge0 = g[ee0];

            // Loop over branches of outgoing edges.
            BGL_FORALL_OUTEDGES(v, e1, g, G) {

                // Reverse complement e1 if necessary.
                edge_descriptor ee1 = e1;
                bool reverse1 = false;
                if(idIsGreaterThanReverseComplement(e1)) {
                    ee1 = g[e1].reverseComplement;
                    reverse1 = true;
                }
                const E& edge1 = g[ee1];

                // Avoid writing links twice.
                if(g[ee0].id > g[ee1].id) {
                    continue;
                }
                if(ee0 == ee1 && reverse0) {
                    continue;
                }

                for(uint64_t i0=0; i0<edge0.ploidy(); i0++) {
                    for(uint64_t i1=0; i1<edge1.ploidy(); i1++) {

                        // To make Bandage happy, we write a Cigar string
                        // consisting of 0M rather than an empty Cigar string.
                        gfa << "L\t" <<
                            edge0.pathId(i0) << "\t" << (reverse0 ? "-" : "+") << "\t" <<
                            edge1.pathId(i1) << "\t" << (reverse1 ? "-" : "+") << "\t0M\n";

                    }
                }
            }
        }
    }
}



string AssemblyGraph2Edge::colorByPeriod(uint64_t branchId) const
{
    if(isBubble()) {
        if(period == 0) {
            return "Green";
        }
        if(branchId == strongestBranchId) {
            return "Grey";
        }
        switch(period) {
        case 2:
            return "Red";
        case 3:
            return "Orange";
        case 4:
            return "Purple";
        default:
            return "Brown";
        }
    } else {
        return "Grey";
    }

}



// For each edge, compute the number of raw sequence bases
// transfered in each direction for gfa output.
void AssemblyGraph2::countTransferredBases()
{
    G& g = *this;

    BGL_FORALL_EDGES(e, g, G) {
        E& edge = g[e];
        edge.backwardTransferCount = 0;
        edge.forwardTransferCount = 0;

        // In this loop, only to it for the edges with it
        // not geater than their reverse complement.
        if(idIsGreaterThanReverseComplement(e)) {
            continue;
        }

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



    // For edges with id greater than their reverse complement, get it
    // from the reverse complemented edge.
    BGL_FORALL_EDGES(e, g, G) {

        if(not idIsGreaterThanReverseComplement(e)) {
            continue;
        }

        E& edge = g[e];
        const edge_descriptor eRc = edge.reverseComplement;
        const E& edgeRc = g[eRc];

        edge.backwardTransferCount = edgeRc.forwardTransferCount;
        edge.forwardTransferCount = edgeRc.backwardTransferCount;
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

        /*const BubbleGraph::vertex_descriptor v = */ add_vertex(BubbleGraphVertex(e, edge), bubbleGraph);
    }

    bubbleGraph.createOrientedReadsTable(readCount);
    bubbleGraph.createEdges();
    cout << "The bubble graph has " << num_vertices(bubbleGraph) <<
        " vertices and " << num_edges(bubbleGraph) << " edges." << endl;

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


    // Remove edges with too few reads.
    const uint64_t minCount = 3;
    vector<BubbleGraph::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        if(bubbleGraph[e].totalCount() < minCount) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const BubbleGraph::edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, bubbleGraph);
    }


    ofstream csv("BubbleGraphEdges.csv");
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

