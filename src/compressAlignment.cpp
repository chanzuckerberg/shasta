#include "compressAlignment.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

// Compress/decompress an alignment to bytes.
// See compressAlignment.hpp for more details.


void shasta::compress(const Alignment&, string&)
{
    SHASTA_ASSERT(0);
}
void shasta::decompress(span<const char>, Alignment)
{
    SHASTA_ASSERT(0);
}

