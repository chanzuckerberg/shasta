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
    uint8_t value;

    // Constructors.
    // To avoid confusion, use a second argument to specify
    // whether the argument is the character representing the base (A, C, G, T),
    // or the corresponding integer value (0 through 3).
    class FromInteger {};
    Base() : value(0) {}
    static Base fromCharacter(char c)
    {
        Base base;
        base.value = BaseInitializer::table[uint8_t(c)];
        if(base.value == 255) {
            throw runtime_error("Invalid base character " + to_string(c));
        }
        return base;
    }
    Base(uint8_t value, FromInteger) : value(value)
    {
    }

    // Return the character representing the base.
    // This always returns an uppercase character,
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
        return Base(uint8_t(3-value), FromInteger());
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
