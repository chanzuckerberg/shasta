// Shasta.
#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



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
