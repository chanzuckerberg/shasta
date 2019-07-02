#include "ShortBaseSequence.hpp"
#include "algorithm.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace ::shasta;

#include <iomanip>



void shasta::testShortBaseSequence()
{
    ShortBaseSequence8 s;
    s.set(0, Base::fromCharacter('T'));
    s.set(1, Base::fromCharacter('C'));
    s.set(2, Base::fromCharacter('G'));
    s.set(3, Base::fromCharacter('T'));
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
    SHASTA_ASSERT(t == s);
}
