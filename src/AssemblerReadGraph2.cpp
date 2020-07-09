// This file contains code for ReadGraph.creationMethod 2.

// Shasta.
#include "Assembler.hpp"
using namespace shasta;


void Assembler::createReadGraph2()
{

    // Mark all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentData.size(), false);

    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);

    SHASTA_ASSERT(0);
}
