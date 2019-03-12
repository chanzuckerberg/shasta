#ifndef CZI_SHASTA_ASSEMBLED_SEGMENT_HPP
#define CZI_SHASTA_ASSEMBLED_SEGMENT_HPP

#include "Base.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class AssembledSegment;
    }
}



// Class to describe a sequence segment assembled
// from an edge of the assemblygraph.
class ChanZuckerberg::shasta::AssembledSegment {
public:

    // The assembled run-length sequence  and repeat counts.
    vector<Base> runLengthSequence;
    vector<uint32_t> repeatCounts;

    // Put back into default-constructed state
    // (except for vector capacities).
    void clear();
};




#endif

