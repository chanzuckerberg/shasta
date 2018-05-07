#include "Base.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;


void ChanZuckerberg::Nanopore2::testBase()
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
