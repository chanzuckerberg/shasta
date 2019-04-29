#include "Assembler.hpp"
#include "computeRunLengthRepresentation.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Analyze assembled sequence by mapping it to
// a given reference sequence used as ground truth.
// The first one is python-callable.
// The edgeId is the same as in FASTA or GFA output
// Begin and end specify the portion of the assembled
// segment to be analyzed and they are expressed in
// raw coordinates (not run-length).
// The reference sequence is given in its raw
// representation (not run-length).
// Assembled sequenceis mapped to reference sequence
// using SeqAn's implementation of the Needleman-Wunsch
// algorithm, which is very expensive in time and memory
// when the sequences are long.

// This is incomplete and unused.

void Assembler::analyzeAssembledSequence(
    AssemblyGraph::EdgeId edgeId,
    uint32_t begin,
    uint32_t end,
    const string& referenceSequenceString
)
{
    const size_t n = referenceSequenceString.size();
    vector<Base> referenceSequence(n);
    for(size_t i=0; i<n; i++) {
        try {
            referenceSequence[i] = Base::fromCharacter(referenceSequenceString[i]);
        } catch(runtime_error&) {
            throw runtime_error("Invalid character " + to_string(int(referenceSequenceString[i])) +
                " at position " + to_string(i) + " in reference sequence.");
        }
    }
    analyzeAssembledSequence(edgeId, begin, end, referenceSequence);
}



void Assembler::analyzeAssembledSequence(
    AssemblyGraph::EdgeId edgeId,
    uint32_t begin,
    uint32_t end,
    const vector<Base>& rawReferenceSequence
)
{
    // Convert the reference sequence to run-length representation.
    vector<Base> runLengthReferenceSequence;
    vector<uint8_t> referenceRepeatCount;
    CZI_ASSERT(
    computeRunLengthRepresentation(
        rawReferenceSequence,
        runLengthReferenceSequence,
        referenceRepeatCount));
}
