#include "Base.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

#include "algorithm.hpp"



BaseInitializer BaseInitializer::singleton;
std::array<uint8_t, 256> BaseInitializer::table;
BaseInitializer::BaseInitializer()
{
    fill(table.begin(), table.end(), 255);
    table['A'] = 0;
    table['C'] = 1;
    table['G'] = 2;
    table['T'] = 3;
    table['a'] = 0;
    table['c'] = 1;
    table['g'] = 2;
    table['t'] = 3;
}


void ChanZuckerberg::shasta::testBase()
{
    const Base A('A', Base::FromCharacter());
    if(A.value) throw runtime_error("A is not 0.");

    const Base C('C', Base::FromCharacter());
    if(C.value != 1) throw runtime_error("C is not 1.");

    const Base G('G', Base::FromCharacter());
    if(G.value != 2) throw runtime_error("G is not 2.");

    const Base T('T', Base::FromCharacter());
    if(T.value != 3) throw runtime_error("T is not 3.");

    cout << A << C << G << T << endl;
}
