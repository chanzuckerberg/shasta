#include "Base.hpp"
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



AlignedBaseInitializer AlignedBaseInitializer::singleton;
std::array<uint8_t, 256> AlignedBaseInitializer::table;
AlignedBaseInitializer::AlignedBaseInitializer()
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
    table['-'] = 4;
}



void shasta::testBase()
{
    const Base A = Base::fromCharacter('A');
    if(A.value != 0) throw runtime_error("A is not 0.");

    const Base C = Base::fromCharacter('C');
    if(C.value != 1) throw runtime_error("C is not 0.");

    const Base G = Base::fromCharacter('G');
    if(G.value != 2) throw runtime_error("G is not 0.");

    const Base T = Base::fromCharacter('T');
    if(T.value != 3) throw runtime_error("R is not 0.");

    cout << A << C << G << T << endl;
}
