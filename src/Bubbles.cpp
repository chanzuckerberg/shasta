#include "Bubbles.hpp"
#include "Assembler.hpp"
#include "prefixLength.hpp"
using namespace shasta;


Bubbles::Bubbles(
    const Assembler& assembler) :
    assembler(assembler)
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
        bubbles.push_back(Bubble());

    }

    cout << suppressedDueToRepeats << " bubbles were suppressed due to 2- or 3-base repeats." << endl;
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

