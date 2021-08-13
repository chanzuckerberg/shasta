#ifndef SHASTA_ALIGNMENT_CANDIDATES_HPP
#define SHASTA_ALIGNMENT_CANDIDATES_HPP

// Shasta.
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"

namespace shasta {
    class AlignmentCandidates;
    class OrientedReadPair;
}



// Alignment candidates found by the LowHash algorithm.
// They all have readId0<readId1.
// They are interpreted with readId0 on strand 0.
class shasta::AlignmentCandidates {
public:
    MemoryMapped::Vector<OrientedReadPair> candidates;

    // For each alignment candidate, we also store a vector of
    // pairs (ordinal0, ordinal1), each containing
    // ordinals in the two oriented reads where identical
    // features (sequences of m markers) were found by the LowHash algorithm.
    // This is only created when using findAlignmentCandidatesLowHashNew
    // (class LowHashNew).
    // This has a vector for each entry in the candidates vector above
    // and is indexed in the same way.
    MemoryMapped::VectorOfVectors< array<uint32_t, 2>, uint64_t> featureOrdinals;

    // The candidate table stores the read pair that each oriented read is involved in.
    // Stores, for each OrientedReadId, a vector of indexes into the alignmentCandidate vector.
    // Indexed by OrientedReadId::getValue(),
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> candidateTable;

    // Method to perform the indexing that fills candidateTable. Allocation of the memory mapped vector requires
    // knowing the number of reads participating in the candidate pairs, and a name + page size to initialize with.
    void computeCandidateTable(ReadId readCount, string largeDataName, size_t largeDataPageSize);

    void unreserve() {
        candidates.unreserve();
        // featureOrdinals is not used by LowHash0
        if (featureOrdinals.isOpenWithWriteAccess()) featureOrdinals.unreserve();
    }

    void clear() {
        candidates.clear();
        // featureOrdinals is not used by LowHash0
        if (featureOrdinals.isOpenWithWriteAccess()) featureOrdinals.clear();
        unreserve();
    }
};

#endif

