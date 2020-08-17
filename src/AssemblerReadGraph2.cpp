#include "Assembler.hpp"
using namespace shasta;



void Assembler::createReadGraph2(
    uint32_t maxAlignmentCount,
    uint32_t maxTrim)
{
    // This boilerplate code creates the read graph using
    // all available alignments.
    vector<bool> keepAlignment(alignmentData.size(), true);
    createReadGraphUsingSelectedAlignments(keepAlignment);
}
