#include "Assembler.hpp"
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



    // Use the Bubbles to flag the alignments we want to keep.
    vector<bool> keepAlignment(alignmentData.size(), false);
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const OrientedReadPair orientedReadPair = alignmentData[alignmentId];
        keepAlignment[alignmentId] = bubbles->allowAlignment(orientedReadPair, useClustering);
    }

    const uint64_t allowedAlignmentCount =
        count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << allowedAlignmentCount << " alignments allowed by bubble analysis out of " <<
        alignmentData.size() << " total." << endl;



    // Keep at most maxAlignmentCount alignments for each read.



    // Missing code.
    SHASTA_ASSERT(0);

    cout << timestamp << "createReadGraphMode1 ends." << endl;
}
