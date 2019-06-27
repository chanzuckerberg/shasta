#include "computeRunLengthRepresentation.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Given the raw representation of a sequence, compute its
// run-length representation.
// This returns false if the sequence contains a homopolymer run
// of more than 255 bases, which cannot be represented
// with a one-byte repeat count.
bool ChanZuckerberg::shasta::computeRunLengthRepresentation(
    const vector<Base>& sequence,
    vector<Base>& runLengthSequence,
    vector<uint8_t>& repeatCount)
{
    runLengthSequence.clear();
    repeatCount.clear();

    for(auto it=sequence.begin(); it!=sequence.end(); ) {
        const Base base = *it;
        uint32_t count = 0;
        while(it!=sequence.end() && *it==base) {
            ++it;
            ++count;
            if(count == 256) {
                return false;
            }
        }
        SHASTA_ASSERT(count > 0);
        SHASTA_ASSERT(count <= 255);
        runLengthSequence.push_back(base);
        repeatCount.push_back(uint8_t(count));
    }

    SHASTA_ASSERT(runLengthSequence.size() == runLengthSequence.size());
    return true;

}
