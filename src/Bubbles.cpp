#include "Bubbles.hpp"
#include "Assembler.hpp"
#include "prefixLength.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>


Bubbles::Bubbles(
    const Assembler& assembler) :
    assembler(assembler)
{
    findBubbles();
    writeBubbles();
    cout << "Found " << bubbles.size() << " bubbles." << endl;

    fillOrientedReadsTable();
    writeOrientedReadsTable();

    createBubbleGraph();
    bubbleGraph.writeGraphviz();
    cout << "The bubble graph has " << num_vertices(bubbleGraph) <<
        " vertices and " << num_edges(bubbleGraph) << " edges." << endl;
}





// For now we only consider diploid bubbles, defined using the
// following strict criteria:
// - Source vertex v0 has out-degree 2.
// - Target vertex v1 has in-degree 2.
// - There are two parallel edges eA and eB, both v0->v1.

void Bubbles::findBubbles()
{
    const bool debug = false;

    // Loop over possible choices for v0 (source vertex of the bubble).
    const AssemblyGraph& assemblyGraph = *assembler.assemblyGraphPointer;
    const AssemblyGraph::VertexId vertexCount = assemblyGraph.vertices.size();
    uint64_t suppressedDueToRepeats = 0;
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
        const auto inEdges1 = assemblyGraph.edgesByTarget[v1];
        if(inEdges1.size() != 2) {
            continue;
        }

        // Assemble the RLE sequence for the two edges.
        vector<Base> sequenceA;
        vector<Base> sequenceB;
        assembler.assembleAssemblyGraphEdgeRleStrict(eA, sequenceA);
        assembler.assembleAssemblyGraphEdgeRleStrict(eB, sequenceB);


        // If the RLE sequences are identical, skip it.
        // This can happen due to missing alignments.
        if(sequenceA == sequenceB) {
            continue;
        }

        if(debug) {
            const MarkerGraph::VertexId u0 = assemblyGraph.vertices[v0];
            const MarkerGraph::VertexId u1 = assemblyGraph.vertices[v1];
            cout << "Bubble at marker graph " << u0 << "->" << u1 << "\n";
            copy(sequenceA.begin(), sequenceA.end(), ostream_iterator<Base>(cout));
            cout << "\n";
            copy(sequenceB.begin(), sequenceB.end(), ostream_iterator<Base>(cout));
            cout << "\n";
        }

        // If the two sequences differ by a 2- of 3-base repeat,
        // this is probably an erroneous bubble.
        if(isShortRepeatCopyNumberDifference(sequenceA, sequenceB)) {
            ++suppressedDueToRepeats;
            if(debug) {
                cout << "Is short repeat difference - skipping.\n";
            }
            continue;
        }

        // We have a diploid bubble.
        bubbles.push_back(Bubble(v0, v1, eA, eB, assembler.markerGraph, assemblyGraph));

    }

    cout << suppressedDueToRepeats << " bubbles were suppressed due to 2- or 3-base repeats." << endl;
}



void Bubbles::writeBubbles()
{
    ofstream csv("Bubbles.csv");
    csv << "BubbleId,VertexId0,VertexId1,CoverageA,CoverageB\n";

    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        const Bubble& bubble = bubbles[bubbleId];
        csv << bubbleId << ",";
        csv << bubble.mv0 << ",";
        csv << bubble.mv1 << ",";
        csv << bubble.orientedReadIds[0].size() << ",";
        csv << bubble.orientedReadIds[1].size() << "\n";
    }
}


// Figure out if two sequences differ only by copy numbers in
// a 2- or 3-base repeat.
bool Bubbles::isShortRepeatCopyNumberDifference(
    const vector<Base>& x,
    const vector<Base>& y)
{
    const bool debug = false;

    // Get the lengths and their differences.
    const uint64_t nx = x.size();
    const uint64_t ny = y.size();

    // When they have the same length, we can certainly return false.
    if(nx == ny) {
        return false;
    }

    // Recursive call so x is shorter than y.
    if(ny < nx) {
        return isShortRepeatCopyNumberDifference(y, x);
    }
    SHASTA_ASSERT(nx < ny);

    // If the length difference is not a multiple of 2 or 3,
    // we can certainly return false.
    const uint64_t dn = ny - nx;
    if( ((dn%2) != 0) and
        ((dn%3) != 0) ) {
        if(debug) {
            cout << "Length difference is not a multiple of 2 or 3\n";
            cout << "nx " << nx << ", ny " << ny << ", dn " << dn << "\n";
        }
        return false;
    }

    const uint64_t commonPrefixLength = shasta::commonPrefixLength(x, y);
    const uint64_t commonSuffixLength = shasta::commonSuffixLength(x, y);

    // Find the portion of y that is not in x.
    uint64_t ix = commonPrefixLength;
    uint64_t iy = commonPrefixLength;
    uint64_t jx = nx - commonSuffixLength;
    uint64_t jy = ny - commonSuffixLength;
    while((jx<ix) or (jy<iy)) {
        ++jx;
        ++jy;
    }

    if(debug) {
        cout << "commonPrefixLength " << commonPrefixLength << "\n";
        cout << "commonSuffixLength " << commonSuffixLength << "\n";
        cout << "nx " << nx << "\n";
        cout << "ny " << ny << "\n";
        cout << "dn " << dn << "\n";
        cout << "ix " << ix << "\n";
        cout << "jx " << jx << "\n";
        cout << "iy " << iy << "\n";
        cout << "jy " << jy << "\n";
    }


    if(ix != jx) {
        // There is more than just an insertion.
        return false;
    }


    // x and y differ by an insertion in iy of range [iy, jy).
    SHASTA_ASSERT(ix == jx);
    SHASTA_ASSERT(jy - iy == dn);



    // Check for k base repeat.
    // We kept the entire common prefix, so we can check just to the left of the insertion.
    for(uint64_t k=2; k<=3; k++) {
        if((dn % k) != 0) {
            continue;
        }
        if(debug) {
            cout << "Trying k = " << k << "\n";
        }

        // Check that the inserted bases are a repeat with period k.
        const uint64_t m = dn / k;
        bool repeatViolationFound = false;
        for(uint64_t i=0; i<m; i++) {
            for(uint64_t j=0; j<k; j++) {
                if(y[iy + i*k + j] != y[iy + j]) {
                    repeatViolationFound = true;
                    break;
                }
            }
        }
        if(repeatViolationFound) {
            // The inserted portion is not a repeat of period k.
            continue;
        }

        // Check the previous k bases in both x and y.
        if(ix < k) {
            continue;
        }
        if(iy < k) {
            continue;
        }
        for(uint64_t j=0; j<k; j++) {
            if(y[iy - k + j] != y[ix + j]) {
                repeatViolationFound = true;
                break;
            }
            if(x[ix - k + j] != y[ix + j]) {
                repeatViolationFound = true;
                break;
            }
        }
        if(repeatViolationFound) {
            continue;
        }

        // It is an insertion in y of a repeat of length k.
        return true;
    }



    // None of the k values we tried worked.
    return false;
}



Bubbles::Bubble::Bubble(
    AssemblyGraph::VertexId av0,
    AssemblyGraph::VertexId av1,
    AssemblyGraph::EdgeId eA,
    AssemblyGraph::EdgeId eB,
    const MarkerGraph& markerGraph,
    const AssemblyGraph& assemblyGraph) :
    av0(av0),
    av1(av1),
    aEdgeIds({eA, eB}),
    mv0(assemblyGraph.vertices[av0]),
    mv1(assemblyGraph.vertices[av1])
{
    fillInOrientedReadIds(markerGraph, assemblyGraph);
}



void Bubbles::Bubble::fillInOrientedReadIds(
    const MarkerGraph& markerGraph,
    const AssemblyGraph& assemblyGraph)
{

    // Loop over both sides of the bubble.
    array<std::set<OrientedReadId>, 2> orientedReadIdSets;
    for(uint64_t i=0; i<2; i++) {
        const AssemblyGraph::EdgeId aEdgeId = aEdgeIds[i];

        // Loop over marker graph edges of this side.
        const span<const MarkerGraph::EdgeId> mEdgeIds = assemblyGraph.edgeLists[aEdgeId];
        for(const MarkerGraph::EdgeId mEdgeId : mEdgeIds) {

            // Loop over MarkerIntervals of this marker graph edge.
            const span<const MarkerInterval> markerIntervals =
                markerGraph.edgeMarkerIntervals[mEdgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                orientedReadIdSets[i].insert(markerInterval.orientedReadId);
            }
        }
    }

    // Now store them, discarding the ones that appear on both sides.
    for(uint64_t i=0; i<2; i++) {
        const std::set<OrientedReadId>& x = orientedReadIdSets[i];
        const std::set<OrientedReadId>& y = orientedReadIdSets[1-i];
        for(const OrientedReadId orientedReadId: x) {
            if(y.find(orientedReadId) == y.end()) {
                orientedReadIds[i].push_back(orientedReadId);
            }
        }
    }
}



void Bubbles::fillOrientedReadsTable()
{
    orientedReadsTable.resize(2* assembler.getReads().readCount());
    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        const Bubble& bubble = bubbles[bubbleId];
        for(uint64_t side=0; side<2; side++) {
            for(const OrientedReadId orientedReadId: bubble.orientedReadIds[side]) {
                orientedReadsTable[orientedReadId.getValue()].push_back(make_pair(bubbleId, side));
            }
        }
    }
}



void Bubbles::writeOrientedReadsTable()
{
    ofstream csv("OrientedReadsTable.csv");
    csv << "OrientedReadId,Bubble,Side\n";
    for(ReadId readId=0; readId<assembler.getReads().readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            for(const auto& v: orientedReadsTable[orientedReadId.getValue()]) {
                csv << orientedReadId << ",";
                csv << v.first << ",";
                csv << v.second << "\n";
            }
        }
    }
}



void Bubbles::createBubbleGraph()
{
    // Create the vertices.
    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        add_vertex(bubbleGraph);
    }


    // Create the edges.
    for(uint64_t bubbleIdA=0; bubbleIdA<bubbles.size(); bubbleIdA++) {
        const Bubble& bubbleA = bubbles[bubbleIdA];

        for(uint64_t sideA=0; sideA<2; sideA++) {
            const vector<OrientedReadId>& orientedReadIds = bubbleA.orientedReadIds[sideA];

            for(const OrientedReadId orientedReadId: orientedReadIds) {
                const auto& v = orientedReadsTable[orientedReadId.getValue()];

                for(const auto& p: v) {
                    const uint64_t bubbleIdB = p.first;
                    if(bubbleIdB <= bubbleIdA) {
                        continue;
                    }
                    const uint64_t sideB = p.second;

                    // Locate the edge between these two bubbles,
                    // creating it if necessary.
                    BubbleGraph::edge_descriptor e;
                    bool edgeWasFound = false;
                    tie(e, edgeWasFound) = edge(bubbleIdA, bubbleIdB, bubbleGraph);
                    if(not edgeWasFound) {
                        tie(e, edgeWasFound) = add_edge(bubbleIdA, bubbleIdB, bubbleGraph);
                    }
                    SHASTA_ASSERT(edgeWasFound);

                    // Update the matrix.
                    BubbleGraphEdge& edge = bubbleGraph[e];
                    ++edge.matrix[sideA][sideB];
                }
            }
        }
    }



    ofstream csv("BubbleGraphEdges.csv");
    csv << "BubbleIdA,BubbleIdB,m00,m11,m01,m10\n";
    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        const BubbleGraphEdge& edge = bubbleGraph[e];
        csv << source(e, bubbleGraph) << ",";
        csv << target(e, bubbleGraph) << ",";
        csv << edge.matrix[0][0] << ",";
        csv << edge.matrix[1][1] << ",";
        csv << edge.matrix[0][1] << ",";
        csv << edge.matrix[1][0] << "\n";
    }
}



void Bubbles::BubbleGraph::writeGraphviz()
{
    const BubbleGraph& bubbleGraph = *this;

    ofstream out("BubbleGraph.dot");
    out << "graph G{\n"
        "node [shape=point];\n";

    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        const BubbleGraphEdge& edge = bubbleGraph[e];
        const uint64_t bubbleIdA = source(e, bubbleGraph);
        const uint64_t bubbleIdB = target(e, bubbleGraph);

        const uint64_t diagonal = edge.matrix[0][0] + edge.matrix[1][1];
        const uint64_t offDiagonal = edge.matrix[0][1] + edge.matrix[1][0];
        const uint64_t total = diagonal + offDiagonal;

        // The concordantRatio is 1 if all is good and 0.5 in the worst case.
        const double concordantRatio = double(max(diagonal, offDiagonal)) / double(total);
        const double hue = (concordantRatio - 0.5) * 0.6666667;

        out << bubbleIdA << "--" << bubbleIdB <<
            " ["
            // "penwidth=" << 0.1*double(total) <<
            " color=\"" << hue << ",1.,0.8\""
            "];\n";
    }

    out << "}\n";
}
