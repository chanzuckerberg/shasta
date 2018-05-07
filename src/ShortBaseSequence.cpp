#include "ShortBaseSequence.hpp"
#include "algorithm.hpp"
#include "CZI_ASSERT.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

#include <iomanip>



void ChanZuckerberg::Nanopore2::testShortBaseSequence()
{
    ShortBaseSequence8 s;
    s.set(0, Base('T', Base::FromCharacter()));
    s.set(1, Base('C', Base::FromCharacter()));
    s.set(2, Base('G', Base::FromCharacter()));
    s.set(3, Base('T', Base::FromCharacter()));
    cout << s << " " << s.reverseComplement(6) << endl;
    s.shiftLeft();
    cout << s << endl;

    // const auto oldFill = cout.fill('0');
    for(const uint8_t x: s.data) {
        cout << std::setw(2) << std::hex << int(x) << endl;
        // cout << int(x) << endl;
    }
    // cout.fill(oldFill);

    // Check that constructor from id does the inverse of function id().
    const ShortBaseSequence8 t(s.id(4), 4);
    CZI_ASSERT(t == s);
}
