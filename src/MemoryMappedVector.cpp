#include "MemoryMappedVector.hpp"


void ChanZuckerberg::shasta::testMemoryMappedVector()
{
    MemoryMapped::Vector<int> x;
    x.createNew("abc", 5);
    x[4] = 18;
    CZI_ASSERT(x[4] == 18);
    x.resize(2*1024*1024 + 20);
    x[2*1024*1024 + 4] = 19;
    CZI_ASSERT(x[2*1024*1024 + 4] == 19);
    x.remove();
}
