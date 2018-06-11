#ifndef CZI_NANOPORE2_MARKER_FINDER_HPP
#define CZI_NANOPORE2_MARKER_FINDER_HPP

// Nanopore2
#include "Marker.hpp"
#include "MultitreadedObject.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {
        class MarkerFinder;
        class LongBaseSequences;
        namespace MemoryMapped {
            template<class T> class Vector;
            template<class Int, class T> class VectorOfVectors;
        }
    }
}



class ChanZuckerberg::Nanopore2::MarkerFinder :
    public MultithreadedObject<MarkerFinder>{
public:

    // The constructor does all the work.
    MarkerFinder(
        size_t k,
        const MemoryMapped::Vector<KmerInfo>& kmerTable,
        LongBaseSequences& reads,
        MemoryMapped::VectorOfVectors<CompressedMarker0, uint64_t>& markers0,
        MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        size_t threadCount);

private:

    // The arguments passed to the constructor.
    size_t k;
    const MemoryMapped::Vector<KmerInfo>& kmerTable;
    LongBaseSequences& reads;
    MemoryMapped::VectorOfVectors<CompressedMarker0, uint64_t>& markers0;
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    size_t threadCount;

    void threadFunction(size_t threadId);

    // In pass 1, we count the number of markers for each
    // read and call reads.incrementCountMultithreaded.
    // In pass 2, we store the markers.
    size_t pass;

};

#endif
