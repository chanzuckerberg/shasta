#ifndef SHASTA_MARKER_FINDER_HPP
#define SHASTA_MARKER_FINDER_HPP

#include "Marker.hpp"
#include "MultithreadedObject.hpp"
#include "Reads.hpp"

namespace shasta {
    class MarkerFinder;
    class LongBaseSequences;

    namespace MemoryMapped {
        template<class T> class Vector;
        template<class Int, class T> class VectorOfVectors;
    }
}



class shasta::MarkerFinder :
    public MultithreadedObject<MarkerFinder>{
public:

    // The constructor does all the work.
    MarkerFinder(
        size_t k,
        const MemoryMapped::Vector<KmerInfo>& kmerTable,
        const Reads& reads,
        MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        size_t threadCount);

private:

    // The arguments passed to the constructor.
    size_t k;
    const MemoryMapped::Vector<KmerInfo>& kmerTable;
    const Reads& reads;
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    size_t threadCount;

    void threadFunction(size_t threadId);

    // In pass 1, we count the number of markers for each
    // read and call reads.incrementCountMultithreaded.
    // In pass 2, we store the markers.
    size_t pass;

};

#endif
