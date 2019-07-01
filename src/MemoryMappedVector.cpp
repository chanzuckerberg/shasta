#include "MemoryMappedObject.hpp"
#include "MemoryMappedVector.hpp"

namespace shasta {
    class MemoryMappedObjectTest {
    public:
        int a;
        int b;
    };
}



void shasta::testMemoryMappedVector()
{
#if 0
    MemoryMapped::Vector<int> x;
    x.createNew("", 2*1024*1024, 5);
    x[4] = 18;
    SHASTA_ASSERT(x[4] == 18);
    x.resize(2*1024*1024 + 20);
    x[2*1024*1024 + 4] = 19;
    SHASTA_ASSERT(x[2*1024*1024 + 4] == 19);
    x.remove();
#endif

#if 1
    MemoryMapped::Object<MemoryMappedObjectTest> x;
    x.createNew("", 2*1024*1024);
    x->a = 2;
    x->b = 3;
    SHASTA_ASSERT(x->a == 2);
    SHASTA_ASSERT(x->b == 3);
#endif
}
