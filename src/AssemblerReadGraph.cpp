/*******************************************************************************

The read graph is a bidirected graph in which each vertex
represents a read.

Using a bidirected graph, we only need
one vertex per read, not two, like in the naive approach,
which uses one vertex for each of the two orientations
of a read.

For more information on the bidirected approach, see
Medvedev, Georgiou, Myers, and Brud,
"Computability of Models for Sequence Assembly",
International Workshop on Algorithms in Bioinformatics, pp 289-301 (2007).
http://medvedevgroup.com/papers/wabi07.pdf
In particular, see Fig. 1 and Section 2.2.

In summary:
- Each edge has two directions ("arrows"), one for each end.
- In a valid path, the two edges into and from a given vertex
  must have opposite orientations at that vertex.
- If the two edge orientations out of a vertex agree
  with the direction of the path,
  the read associated with the vertex is on strand 0 (unchanged).
  Otherwise, it is on strand 1 (reverse complemented).

We store edges with the lowest numbered read as the first read.

We use the standard approach to construct string graphs
(Myers, "The fragment assembly string graph" (2005),
doi:10.1093/bioinformatics/bti111,
http://www.cs.utoronto.ca/~brudno/csc2427/myers.pdf),
which prescribes that contained reads should not be included in the graph.
A contained read is a read that has an alignment
with a longer read covering the entire read.

However, in order not to lose coverage from contained reads
(which would result in a unacceptable reduction in coverage),
for each contained read we keep track of the containing read
that achieves the best alignment. This information is used later,
to incorporate contained reads in the marker graph.

*******************************************************************************/


// Shasta.
#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace shasta;


#if 0
// Create the global read graph.
void Assembler::createReadGraph(uint32_t maxTrim)
{
    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    CZI_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;



    // Find which oriented reads are contained.
    vector<bool> isContained(orientedReadCount);
    uint32_t containedCount = 0;
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {

            // Get the information we need, including the alignments for this oriented read.
            const OrientedReadId orientedReadId0(readId0, strand0);
            const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
                findOrientedAlignments(orientedReadId0);
            const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());

            // Check all the alignments.
            isContained[orientedReadId0.getValue()] = false;
            for(const auto& p: alignments) {
                const AlignmentInfo& alignmentInfo = p.second;
                const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
                const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
                if(leftTrim0<=maxTrim && rightTrim0<=maxTrim) {
                    isContained[orientedReadId0.getValue()] = true;
                    ++containedCount;
                    break;
                }
            }
        }
    }

    cout << "Found " << containedCount << " contained oriented reads out of ";
    cout << orientedReadCount << " total." << endl;
    cout << "Number of non-contained oriented reads is ";
    cout << orientedReadCount - containedCount << "." << endl;


    ofstream graphOut("ReadGraph.dot");
    graphOut << "graph G {\n";
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);
            if(isContained[orientedReadId0.getValue()]) {
                continue;
            }

            // Get the information we need, including the alignments for this oriented read.
            const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
                findOrientedAlignments(orientedReadId0);
            const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());

            // Check all the alignments.
            for(const auto& p: alignments) {
                const OrientedReadId orientedReadId1 = p.first;
                if(isContained[orientedReadId1.getValue()]) {
                    continue;
                }
                if(orientedReadId1 < orientedReadId0) {
                    continue;
                }

                const uint32_t markerCount1 = uint32_t(markers[orientedReadId0.getValue()].size());
                const AlignmentInfo& alignmentInfo = p.second;

                // If there is too much trim on both sides of orientedReadId0,
                // this is is a suspicious alignment. Discard.
                const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
                const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
                if(leftTrim0 > maxTrim && rightTrim0 > maxTrim) {
                    continue;
                }

                // If there is too much trim on both sides of orientedReadId1,
                // this is is a suspicious alignment. Discard.
                const uint32_t leftTrim1 = alignmentInfo.firstOrdinals.second;
                const uint32_t rightTrim1 = markerCount1 - 1 - alignmentInfo.lastOrdinals.second;
                if(leftTrim1 > maxTrim && rightTrim1 > maxTrim) {
                    continue;
                }

                graphOut << orientedReadId0.getValue() << "--";
                graphOut << orientedReadId1.getValue() << ";\n";
            }
        }
    }
    graphOut << "}\n";
}
#else



// Create the global read graph.
void Assembler::createReadGraph(uint32_t maxTrim)
{
    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    CZI_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Allocate containingOrientedReadId and initialize it to
    // OrientedReadId::invalid().
    // This marks all read as not contained.
    containingOrientedReadId.createNew(largeDataName("ContainingOrientedReadId"), largeDataPageSize);
    containingOrientedReadId.resize(readCount);
    fill(containingOrientedReadId.begin(), containingOrientedReadId.end(), OrientedReadId::invalid());


    // Find which reads are contained.
    // For the ones that are contained, set containingOrientedReadId
    // to OrientedReadId(0, 0), for now.
    uint32_t containedCount = 0;
    for(ReadId readId0=0; readId0<readCount; readId0++) {

        // Get the information we need.
        const Strand strand0 = 0;
        const OrientedReadId orientedReadId0(readId0, strand0);
        const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
            findOrientedAlignments(orientedReadId0);
        const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());

        // Check all the alignments.
        for(const auto& p: alignments) {
            const AlignmentInfo& alignmentInfo = p.second;
            const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
            const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
            if(leftTrim0<=maxTrim && rightTrim0<=maxTrim) {
                containingOrientedReadId[readId0] = OrientedReadId(0, 0);
                ++containedCount;
                break;
            }
        }
    }
    cout << "Found " << containedCount << " contained reads out of ";
    cout << readCount << " total." << endl;
    cout << "Number of non-contained reads is ";
    cout << readCount - containedCount << "." << endl;



    // For each contained read, find the best containing oriented read.
    // This is read that achieves the best alignment -
    // that is, the alignment with the greatest number of markers.
    // Note that this best containing oriented read could
    // also itself be contained. This is possible because
    // we are not guaranteed to have all the alignments.
    // That is, if B contains A and C contains B,
    // we don't necessarily have the alignment between A and C.
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        if(!isContainedRead(readId0)) {
            continue;
        }

        // Get the information we need.
        const Strand strand0 = 0;
        const OrientedReadId orientedReadId0(readId0, strand0);
        const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
            findOrientedAlignments(orientedReadId0);
        const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());

        // Check all the alignments.
        OrientedReadId bestContaining = OrientedReadId::invalid();
        uint32_t bestMarkerCount = 0;
        for(const auto& p: alignments) {
            const OrientedReadId orientedReadId1 = p.first;
            const AlignmentInfo& alignmentInfo = p.second;
            const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
            const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
            if(leftTrim0<=maxTrim && rightTrim0<=maxTrim) {
                if(alignmentInfo.markerCount > bestMarkerCount) {
                    bestContaining = orientedReadId1;
                    bestMarkerCount = alignmentInfo.markerCount;
                }
            }
        }
        CZI_ASSERT(bestContaining != OrientedReadId::invalid());
        containingOrientedReadId[readId0] = bestContaining;
    }
}

#endif



void Assembler::accessReadGraph()
{
    containingOrientedReadId.accessExistingReadOnly(largeDataName("ContainingOrientedReadId"));
}



// Follow the chain of containing reads until we reach a non-contained read.
OrientedReadId Assembler::findContainingReadRecursive(OrientedReadId orientedReadId) const
{
    while(true) {
        const ReadId readId = orientedReadId.getReadId();
        if(!isContainedRead(readId)) {
            return orientedReadId;
        }
        OrientedReadId containing = containingOrientedReadId[readId];
        if(orientedReadId.getStrand() == 1) {
            containing.flipStrand();
        }
        orientedReadId = containing;
    }
}


