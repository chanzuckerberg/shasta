#include "testMarginCore.hpp"
#include "CZI_ASSERT.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

#if 0
#define delete deleteNotAKeyword
#include "marginPhase/callConsensus.h"
#undef delete
#endif

#include "stdexcept.hpp"
#include "string.hpp"
#include "vector.hpp"



void ChanZuckerberg::shasta::testMarginCore()
{
#if 0
    // Get the parameters.
    const string fileName = "FileName";
    PolishParams* parameters = getConsensusParameters(
        const_cast<char*>(fileName.c_str()));
    if(!parameters) {
        throw runtime_error("testMarginCore: Error opening " + fileName);
    }

    // Vectors to contain the input to callConsensus.
    vector<string> sequences;
    vector< vector<uint8_t> > repeatCounts;
    static_assert(sizeof(bool) == sizeof(uint8_t), "Unexpected sizeof(bool) in testMarginCore.");
    vector<uint8_t> strands;

    // Fill them in
    sequences.push_back("ACTG");
    repeatCounts.push_back(vector<uint8_t>({1, 1, 1, 1}));
    strands.push_back(uint8_t(false));

    // Get the read count.
    const size_t readCount = sequences.size();

    // Sanity checks.
    CZI_ASSERT(repeatCounts.size() == readCount);
    CZI_ASSERT(strands.size() == readCount);

    // Create vectors of pointers to be passed to callConsensus.
    vector<char*> sequencePointers(readCount);
    vector<uint8_t*> repeatCountPointers(readCount);
    for(size_t i=0; i<readCount; i++) {
        sequencePointers[i] = const_cast<char*>(sequences[i].data());
        repeatCountPointers[i] = const_cast<uint8_t*>(repeatCounts[i].data());
    }

    // Do the work.
    callConsensus(
        readCount,
        sequencePointers.data(),
        repeatCountPointers.data(),
        reinterpret_cast<bool*>(strands.data()),
        parameters);

    // Destroy the parameters.
    destroyConsensusParameters(parameters);
#endif
}
