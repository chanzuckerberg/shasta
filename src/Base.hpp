#ifndef CZI_NANOPORE2_BASE_HPP
#define CZI_NANOPORE2_BASE_HPP


#include "array.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {
        class Base;
        class BaseInitializer;
        inline ostream& operator<<(ostream&, Base);
        void testBase();
    }
}


// Class used only to store a static look up table
// use by the Base constructors toc convert
// chacters to bases.
class ChanZuckerberg::Nanopore2::BaseInitializer{
public:
    BaseInitializer();
    static array<uint8_t, 256> table;
    static BaseInitializer singleton;
};



// Describes a single base.
// Represented as 1=byte integer as: A=0, C=1, G=2, T=3.
// This choice of representation facilitates the computation of the complement.
class ChanZuckerberg::Nanopore2::Base {
public:
    uint8_t value;

    // Constructors.
    // To avoid confusion, use a second argument to specify
    // whether the argument is the character representing the base (A, C, G, T),
    // or the corresponding integer value (0 through 3).
    class FromCharacter {};
    class FromCharacterNoException {};
    class FromInteger {};
    Base() : value(0) {}
    Base(char c, FromCharacter)
    {
        value = BaseInitializer::table[uint8_t(c)];
        if(value == 255) {
            throw runtime_error("Invalid base character " + to_string(c));
        }
    }
    Base(char c, FromCharacterNoException)
    {
        value = BaseInitializer::table[uint8_t(c)];
    }
    Base(uint8_t value, FromInteger) : value(value)
    {
    }

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



inline std::ostream& ChanZuckerberg::Nanopore2::operator<<(
    std::ostream& s,
    ChanZuckerberg::Nanopore2::Base base)
{
    s << base.character();
    return s;
}





#endif
