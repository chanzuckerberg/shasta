#ifndef SHASTA_READ_ID_HPP
#define SHASTA_READ_ID_HPP

// Shasta.
#include "SHASTA_ASSERT.hpp"
#include "shastaTypes.hpp"

// Standard libraries.
#include <cstdlib>
#include "cstdint.hpp"
#include "iostream.hpp"
#include <limits>
#include "stdexcept.hpp"
#include "string.hpp"



namespace shasta {

    const ReadId invalidReadId = std::numeric_limits<ReadId>::max();

    // Class used to identify an oriented read,
    // that is a read, possibly reverse complemented.
    class OrientedReadId;
    inline ostream& operator<<(ostream&, OrientedReadId);

}


// Class used to identify an oriented read,
// that is a read, possibly reverse complemented.
// The strand stored in the least significant
// bit is 0 if the oriented read is identical
// to the original read and 1 if it is reverse complemented
class shasta::OrientedReadId {
public:
    OrientedReadId() : value(std::numeric_limits<ReadId>::max()) {}
    OrientedReadId(ReadId readId, Strand strand) : value((readId<<1) | strand)
    {
        SHASTA_ASSERT(strand < 2);
    }

    // This constructor is confusing so I am removing it
    // and replacing it with fromValue.
    OrientedReadId(ReadId) = delete;
    static OrientedReadId fromValue(ReadId value)
    {
        OrientedReadId orientedReadId;
        orientedReadId.value = value;
        return orientedReadId;
    }



    // Constructor from a string of the form readId-strand.
    // Throws an exception on failure.
    explicit OrientedReadId(const string& s)
    {
        // Look for the dash.
        const auto dashPosition = s.find_first_of('-');
        if(dashPosition == string::npos) {
            throw runtime_error("Invalid oriented read id - dash is missing: " + s);
        }

        // Check that only one character follows the dash.
        if(dashPosition != s.size()-2) {
            throw runtime_error("Invalid oriented read id - strand portion is too long: " + s);
        }

        // Parse the Strand.
        const char strandCharacter = s[dashPosition+1];
        Strand strand;
        switch(strandCharacter) {
        case '0':
            strand = 0;
            break;
        case '1':
            strand = 1;
            break;
        default:
            throw runtime_error("Invalid oriented read id - invalid strand: " + s);
        }

        // Parse the ReadId.
        const string readIdString = s.substr(0, dashPosition);
        for(const char c: readIdString) {
            if(not std::isdigit(c)) {
                throw runtime_error("Invalid oriented read id - invalid character in ReadId: " + s);
            }
        }
        const ReadId readId = std::atoi(readIdString.c_str());

        // All good. Store it.
        value = (readId<<1) | strand;
    }



    ReadId getReadId() const
    {
        return value >> 1;
    }
    Strand getStrand() const
    {
        return value & 1;
    }

    void flipStrand()
    {
        value ^= Int(1);
    }

    // Return the integer value, which can be used as an index intoi Assembler::orientedReadIds.
    ReadId getValue() const
    {
        return value;
    }

    // Return a string representing this OrientedReadId.
    string getString() const
    {
        return to_string(getReadId()) + "-" + to_string(getStrand());
    }

    // Return an invalid OrientedReadId.
    static OrientedReadId invalid()
    {
        return OrientedReadId();
    }

    bool operator==(const OrientedReadId& that) const
    {
        return value == that.value;
    }
    bool operator!=(const OrientedReadId& that) const
    {
        return value != that.value;
    }
    bool operator<(const OrientedReadId& that) const
    {
        return value < that.value;
    }
    bool operator<=(const OrientedReadId& that) const
    {
        return value <= that.value;
    }
    bool operator>(const OrientedReadId& that) const
    {
        return value > that.value;
    }

    using Int = ReadId;
private:
    Int value;
};



inline std::ostream& shasta::operator<<(
    std::ostream& s,
    OrientedReadId orientedReadId)
{
    s << orientedReadId.getString();
    return s;
}


#endif
