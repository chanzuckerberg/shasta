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

        ordinal0 = ordinals[0];
        ordinal1 = ordinals[1];

        // Find the end of the streak.
        uint32_t n = 1;
        for(uint64_t j=i+1; j<alignment.ordinals.size(); j++, n++) {
            const auto& nextOrdinals = alignment.ordinals[j];
            if(nextOrdinals[0] != ordinal0 + 1) {
                break;
            }
            if(nextOrdinals[1] != ordinal1 + 1) {
                break;
            }
            ++ordinal0;
            ++ordinal1;
        }

        i += n;

        // cout << "(Skip0, Skip1, N) = (" << skip0 << ", " << skip1 << ", " << n << ")" << endl;
        // Serialize skipx, skipy, n into the output string
        // using the most compact format possible.
        if(Format0::ok(skip0, skip1, n)) {
            Format0 format(skip0, skip1, n);
            s.append(reinterpret_cast<char*>(&format), sizeof(format));
        } else if(Format1::ok(skip0, skip1, n)) {
            Format1 format(skip0, skip1, n);
            s.append(reinterpret_cast<char*>(&format), sizeof(format));
        } else if(Format2::ok(skip0, skip1, n)) {
            Format2 format(skip0, skip1, n);
            s.append(reinterpret_cast<char*>(&format), sizeof(format));
        } else if(Format3::ok(skip0, skip1, n)) {
            Format3 format(skip0, skip1, n);
            s.append(reinterpret_cast<char*>(&format), sizeof(format));
        } else if(Format4::ok(skip0, skip1, n)) {
            Format4 format(skip0, skip1, n);
            s.append(reinterpret_cast<char*>(&format), sizeof(format));
        } else {
            throw std::runtime_error("Shasta Internal Error: Unable to compress alignment streak: "
                "skip0 " + to_string(skip0) +
                ", skip1 " + to_string(skip1) +
                ", n " + to_string(n));
        }
    }
}


void shasta::decompress(span<const char> s, Alignment& alignment)
{
    uint32_t ordinal0 = 0;
    uint32_t ordinal1 = 0;
    
    size_t pos = 0;
    while (pos < s.size()) {
        uint32_t streakSize;
        int32_t skip0, skip1;

        uint8_t formatIdentifier = compressAlignment::extractFormatIdentifier(s[pos]);
        if (formatIdentifier == Format0::id) {
            // Format0
            Format0 const *f0 = reinterpret_cast<Format0*>(const_cast<char*>(&s[pos]));
            pos += sizeof(Format0);
            skip0 = f0->skip0;
            skip1 = f0->skip1;
            streakSize = f0->n();
            
        } else if (formatIdentifier == Format1::id) {
            // Format1
            Format1 const *f1 = reinterpret_cast<Format1*>(const_cast<char*>(&s[pos]));
            pos += sizeof(Format1);
            skip0 = f1->skip0;
            skip1 = f1->skip1;
            streakSize = f1->n();
            
        } else if (formatIdentifier == Format2::id) {
            // Format2
            Format2 const *f2 = reinterpret_cast<Format2*>(const_cast<char*>(&s[pos]));
            pos += sizeof(Format2);
            skip0 = f2->skip0;
            skip1 = f2->skip1;
            streakSize = f2->n();

        } else if (formatIdentifier == Format3::id) {
            // Format3
            Format3 const *f3 = reinterpret_cast<Format3*>(const_cast<char*>(&s[pos]));
            pos += sizeof(Format3);
            skip0 = f3->skip0;
            skip1 = f3->skip1;
            streakSize = f3->n();

        } else {
            // Format4
            Format4 const *f4 = reinterpret_cast<Format4*>(const_cast<char*>(&s[pos]));
            pos += sizeof(Format4);
            skip0 = f4->skip0;
            skip1 = f4->skip1;
            streakSize = f4->n();
        }

        // cout << "(Skip0, Skip1, N) = (" << skip0 << ", " << skip1 << ", " << streakSize << ")" << endl;
       
        ordinal0 += skip0;
        ordinal1 += skip1;
        for(uint32_t i = 0; i < streakSize; i++) {
            alignment.ordinals.push_back(array<uint32_t, 2>({ordinal0+i, ordinal1+i}));
        }

        ordinal0 += (streakSize-1);
        ordinal1 += (streakSize-1);
    }
}

uint8_t shasta::compressAlignment::extractFormatIdentifier(const char c) {
    if ((c & Format0::idMask) == Format0::id) {
        return Format0::id;
    } else if ((c & Format1::idMask) == Format1::id) {
        return Format1::id;
    } else if ((c & Format2::idMask) == Format2::id) {
        return Format2::id;
    } else if ((c & Format3::idMask) == Format3::id) {
        return Format3::id;
    } else if ((c & Format4::idMask) == Format4::id) {
        return Format4::id;
    } else {
        string errorStr("Shasta Internal Error: Unsupported format in byte - ");
        for (uint8_t i = 0; i < sizeof(c) * 8; i++) {                                                          
            errorStr.append(1, !!((c << i) & 0x80) ? '1': '0');                                                                     
        }
        errorStr.append(".");
        throw std::runtime_error(errorStr);
    }
}

void shasta::testAlignmentCompression() {
    vector< array<uint32_t, 2> > ordinals {
        {300, 200}, // First streak begins here (Format2)
        {301, 201},
        {302, 202},
        {305, 206}, // Second streak begins here (Format1)
        {306, 207},
        {320, 250}, // Third streak begins here (Format2)
        {321, 251},
        {322, 252},
        {323, 253},
        {325, 255}, // Fourth streak begins here (Format0)
        {326, 256},
        {350, 257}, // Fifth streak begins here (Format2)
        {351, 258},
        {352, 259},
        {353, 260},
        {354, 261},
        {1000, 400}, // Sixth streak begins here (Format3)
        {1001, 401},
        {1002, 402},
        {600000, 500000}, // Seventh streak begins here (Format4)
        {600001, 500001},
        {500000, 500005}, // Eighth streak (negative skip0) (Format3)
        {500001, 500007},
        {500002, 500008},
        {500003, 500009},
        {500004, 500010},
        {500005, 500011},
        {500006, 500012},
        {500007, 500013},
        {500008, 500014},
    };
  
    Alignment alignment;
    alignment.ordinals = ordinals;

    // cout << "---Original alignment-----" << endl;
    // for(const auto& x: alignment.ordinals) {
    //     cout << x[0] << ", " << x[1] << endl;
    // }
    // cout << "--------------------------" << endl;

    string compressed;
    shasta::compress(alignment, compressed);

    cout << "Uncompressed size = " << ordinals.size() * sizeof(uint32_t) * 2 << " bytes." << endl;
    cout << "Compressed size = " << compressed.size() << " bytes." << endl;

    Alignment decompressedAlignment;
    span<const char> spanOfBytes(compressed.c_str(), compressed.c_str() + compressed.size());
    shasta::decompress(spanOfBytes, decompressedAlignment);

    // cout << "---Decompressed Alignment -----" << endl;
    // for(const auto& x: decompressedAlignment.ordinals) {
    //     cout << x[0] << ", " << x[1] << endl;
    // }
    // cout << "-------------------------------" << endl;

    SHASTA_ASSERT(alignment.ordinals == decompressedAlignment.ordinals);
}
