#include "MemoryMappedAllocator.hpp"
using namespace shasta;

#include <list>
#include <map>
#include "utility.hpp"

void shasta::MemoryMapped::testMemoryMappedAllocator()
{
    try {
        // Create the low level byte allocator.
        const uint64_t pageSize = 2*1024*1024;
        ByteAllocator byteAllocator(
            "Data/testMemoryMappedAllocator", pageSize, 8*pageSize);

        // Create a vector, a list, and a map, all using the same allocator.

        cout << "vector" << endl;
        Allocator<uint16_t> vectorAllocator(byteAllocator);
        vector<uint16_t, Allocator<uint16_t>> v(2, vectorAllocator);
        v[0] = 10;
        v[1] = 10;
        for(uint16_t i=0; i<10; i++) {
            v.push_back(i);
        }

        cout << "list" << endl;
        Allocator<char> listAllocator(byteAllocator);
        std::list<char, Allocator<char>> l(listAllocator);
        for(char c=0; c<10; c++) {
            l.push_back(c);
        }

        cout << "map" << endl;
        Allocator< pair<const int, uint16_t> > mapAllocator(byteAllocator);
        std::map<int, uint16_t, std::less<int>, Allocator< pair<const int, uint16_t> > >
            m(mapAllocator);
        m[5] = 10;
        m[20] = 2;

        cout << "vector:" << endl;
        for(const auto& x: v) {
            cout << x << endl;
        }

        cout << "list:" << endl;
        for(const auto& x: l) {
            cout << int(x) << endl;
        }

        cout << "map:" << endl;
        for(const auto& p: m) {
            cout << p.first << " " << p.second << endl;
        }


        // This will make it throw.
        v.resize(100000000);

    } catch(BadAllocation& e) {
        cout << "BadAllocation was thrown." << endl;
    }
}
