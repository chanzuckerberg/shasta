#include "compressAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace compressAlignment;

// Compress/decompress an alignment to bytes.
// See compressAlignment.hpp for more details.



void shasta::compress(const Alignment& alignment, string& s)
{
    s.clear();
    uint32_t ordinal0 = 0;
    uint32_t ordinal1 = 0;

    for(uint64_t i=0; i<alignment.ordinals.size(); ) {

        // A new streak begins here.
        const auto& ordinals = alignment.ordinals[i];

        // Compute the skip values for this streak.
        const int32_t skip0 = int32_t(ordinals[0]) - int32_t(ordinal0);
        const int32_t skip1 = int32_t(ordinals[1]) - int32_t(ordinal1);

        // Find the end of the streak.
        uint32_t n = 0;
        for(uint64_t j=i+1; j<alignment.ordinals.size(); j++, n++) {
            ++ordinal0;
            ++ordinal1;
            const auto& nextOrdinals = alignment.ordinals[j];
            if(nextOrdinals[0] != ordinal0) {
                break;
            }
            if(nextOrdinals[1] != ordinal1) {
                break;
            }
        }


        // Serialize skipx, skipy, n into the output string
        // using the most compact format possible.
        if(Format0::ok(skip0, skip1, n)) {
            Format0 format(skip0, skip1, n);
            s.append(reinterpret_cast<char*>(&format), sizeof(format));

        } else if(Format1::ok(skip0, skip1, n)) {
            Format1 format(skip0, skip1, n);
            s.append(reinterpret_cast<char*>(&format), sizeof(format));
        } else {
            throw runtime_error("Unable to compress alignment streak: "
                "skip0 " + to_string(skip0) +
                ", skip1 " + to_string(skip1) +
                ", n " + to_string(n));
        }
    }
}




void shasta::decompress(span<const char>, Alignment)
{
    SHASTA_ASSERT(0);
}

