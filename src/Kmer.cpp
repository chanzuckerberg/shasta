#include "Kmer.hpp"

using namespace shasta;
using namespace kmer;

// 4^k
uint64_t shasta::kmer::totalCount(uint64_t k) {
    return(1ULL << (2ULL*k));
}