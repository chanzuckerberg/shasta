#include "Assembler.hpp"
using namespace shasta;



// Read graph creation for mode 1 assembly.
void Assembler::createReadGraphMode1(uint64_t maxAlignmentCount)
{
    cout << timestamp << "createReadGraphMode1 begins." << endl;

    // Check gthat we have what we need.
    checkAlignmentDataAreOpen();
    SHASTA_ASSERT(bubbles);


    cout << timestamp << "createReadGraphMode1 ends." << endl;
    SHASTA_ASSERT(0);
}
