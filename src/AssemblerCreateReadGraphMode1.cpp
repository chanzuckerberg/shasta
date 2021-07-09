#include "Assembler.hpp"
#include "Bubbles.hpp"
#include "Reads.hpp"
using namespace shasta;



// Read graph creation for mode 1 assembly.
void Assembler::createReadGraphMode1(uint64_t maxAlignmentCount)
{
    // Parameters that should be exposed when the code stabilizes.
    const bool useClustering = true;

    cout << timestamp << "createReadGraphMode1 begins." << endl;

    // Check that we have what we need.
    checkAlignmentDataAreOpen();
    SHASTA_ASSERT(bubbles);



    // Use the Bubbles to flag the allowed alignments.
    vector<bool> alignmentIsAllowed(alignmentData.size(), false);
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const OrientedReadPair orientedReadPair = alignmentData[alignmentId];
        alignmentIsAllowed[alignmentId] = bubbles->allowAlignment(orientedReadPair, useClustering);
    }

    const uint64_t allowedAlignmentCount =
        count(alignmentIsAllowed.begin(), alignmentIsAllowed.end(), true);
    cout << allowedAlignmentCount << " alignments allowed by bubble analysis out of " <<
        alignmentData.size() << " total." << endl;



    // Mark all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentData.size(), false);



    // Vector to keep the alignments for each read,
    // with their number of markers.
    // Contains pairs(marker count, alignment id).
    vector< pair<uint32_t, uint32_t> > readAlignments;



    // Loop over reads.
    for(ReadId readId=0; readId<getReads().readCount(); readId++) {

        // Gather the allowed alignments for this read, each with its number of markers.
        readAlignments.clear();
        for(const uint32_t alignmentId: alignmentTable[OrientedReadId(readId, 0).getValue()]) {
            if(alignmentIsAllowed[alignmentId]) {
                const AlignmentData& alignment = alignmentData[alignmentId];
                readAlignments.push_back(make_pair(alignment.info.markerCount, alignmentId));
            }
        }

        // Keep the best maxAlignmentCount.
        if(readAlignments.size() > maxAlignmentCount) {
            std::nth_element(
                readAlignments.begin(),
                readAlignments.begin() + maxAlignmentCount,
                readAlignments.end(),
                std::greater< pair<uint32_t, uint32_t> >());
            readAlignments.resize(maxAlignmentCount);
        }

        // Mark the surviving alignments as to be kept.
        for(const auto& p: readAlignments) {
            const uint32_t alignmentId = p.second;
            keepAlignment[alignmentId] = true;
        }
    }
    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Initially keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;


    // Add alignments to avoid coverage holes.
    fixCoverageHoles(keepAlignment);
    cout << "After fixing coverage holes, keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    // Create the read graph using the alignments we selected.
    readGraph.remove();
    createReadGraphUsingSelectedAlignments(keepAlignment);

    cout << timestamp << "createReadGraphMode1 ends." << endl;
}
