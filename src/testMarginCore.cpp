#ifndef SHASTA_STATIC_EXECUTABLE

// Shasta.
#include "testMarginCore.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// MarginCore.
#include "marginPhase/callConsensus.h"

// Standard library.
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"
#include "vector.hpp"



void ChanZuckerberg::shasta::testMarginCore()
{
    // Get the parameters.
    const string fileName = "MarginPhase.json";
    PolishParams* parameters = getConsensusParameters(
        const_cast<char*>(fileName.c_str()));
    if(!parameters) {
        throw runtime_error("testMarginCore: Error opening " + fileName);
    }

    // Vectors to contain the input to callConsensus.
    vector<string> sequences;
    vector< vector<uint8_t> > repeatCounts;
    vector<uint8_t> strands;

    // Fill them in
    sequences.push_back("AGTCG");
    repeatCounts.push_back(vector<uint8_t>({1, 1, 1, 1, 3}));
    strands.push_back(0);
    sequences.push_back("AGTG");
    repeatCounts.push_back(vector<uint8_t>({1, 1, 1, 3}));
    strands.push_back(0);
    sequences.push_back("AGTG");
    repeatCounts.push_back(vector<uint8_t>({1, 1, 1, 3}));
    strands.push_back(0);

    // Get the read count.
    const size_t readCount = sequences.size();

    // Sanity checks.
    SHASTA_ASSERT(repeatCounts.size() == readCount);
    SHASTA_ASSERT(strands.size() == readCount);

    // Create vectors of pointers to be passed to callConsensus.
    vector<char*> sequencePointers(readCount);
    vector<uint8_t*> repeatCountPointers(readCount);
    for(size_t i=0; i<readCount; i++) {
        sequencePointers[i] = const_cast<char*>(sequences[i].data());
        repeatCountPointers[i] = const_cast<uint8_t*>(repeatCounts[i].data());
    }

    // Write out the input.
    cout << "Input to callConsensus:" << endl;
    for(size_t i=0; i<readCount; i++) {
        cout << "Read " << i << ": ";
        cout << sequences[i] << " strand " << int(strands[i]) << ", repeat counts ";
        for(const uint8_t r: repeatCounts[i]) {
            cout << " " << int(r);
        }
        cout << endl;
    }

    // Do the work.
    RleString* consensusPointer = callConsensus(
        readCount,
        sequencePointers.data(),
        repeatCountPointers.data(),
        strands.data(),
        parameters);
    const RleString& consensus = *consensusPointer;
    cout << "Output from callConsensus:" << endl;
    for(int64_t i=0; i<consensus.length; i++) {
        cout << consensus.rleString[i];
    }
    cout << ", repeat counts";
    for(int64_t i=0; i<consensus.length; i++) {
        cout << " " << consensus.repeatCounts[i];
    }
    cout << endl;

    // Clean up.
    destroyRleString(consensusPointer);
    destroyConsensusParameters(parameters);
}

#endif

