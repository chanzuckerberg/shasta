#include "AssembledSegment.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



void AssembledSegment::clear()
{
    runLengthSequence.clear();
    repeatCounts.clear();
}
