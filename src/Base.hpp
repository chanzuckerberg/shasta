#ifndef CZI_SHASTA_BASE_HPP
#define CZI_SHASTA_BASE_HPP


#include "array.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class Base;
        class BaseInitializer;
        inline ostream& operator<<(ostream&, Base);
        void testBase();
    }
}


// Class used only to store a static look up table
// use by the Base::fromCharacter to convert
// characters to bases.
class ChanZuckerberg::shasta::BaseInitializer{
public:
    BaseInitializer();
    static array<uint8_t, 256> table;
    static BaseInitializer singleton;
};



// Describes a single base.
// Represented as 1=byte integer as: A=0, C=1, G=2, T=3.
// This choice of representation facilitates the computation of the complement.
class ChanZuckerberg::shasta::Base {
public:

    // Tha byte value is always one of 0, 1, 2, 3.
    uint8_t value;

    // The default constructor constructs A.
    Base() : value(0) {}

    // We use static member functions instead of constructors.
    // This is safer due to the possibility of unwanted
    // conversions between characters and integers,
    // or confusion between the value stored (0, 1, 2, 3) and
    // the representing character (A, C, G, T).

    // Construct from a character.
    static Base fromCharacter(char c)
    {
        Base base;
        base.value = BaseInitializer::table[uint8_t(c)];
        if(base.value == 255) {
            throw runtime_error("Invalid base character " + to_string(c));
        }
        return base;
    }

    // Construct from an integer.
    // This does not check validity.
    static Base fromInteger(uint8_t i)
    {
        Base base;
        base.value = i;
        return base;
    }
    static Base fromInteger(uint16_t i)
    {
        Base base;
        base.value = uint8_t(i);
        return base;
    }
    static Base fromInteger(uint32_t i)
    {
        Base base;
        base.value = uint8_t(i);
        return base;
    }
    static Base fromInteger(uint64_t i)
    {
        Base base;
        base.value = uint8_t(i);
        return base;
    }

    // Return the character representing the base.
    // This always returns an upper case character,
    // regardless of how the base was constructed.
    char character() const
    {
        switch(value) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default:
            throw runtime_error("Invalid base value " + to_string(value));
        }
    }

    // Return the complement of this base, without changing it.
    Base complement() const
    {
        return fromInteger(uint8_t(3 - int(value)));
    }

    // Replace this base with its complement.
    void complementInPlace()
    {
        value = (~value) & 3;
    }

    bool operator==(Base that) const
    {
        return value == that.value;
    }
    bool operator<(Base that) const
    {
        return value < that.value;
    }
};



inline std::ostream& ChanZuckerberg::shasta::operator<<(
    std::ostream& s,
    ChanZuckerberg::shasta::Base base)
{
    s << base.character();
    return s;
}





#endif
