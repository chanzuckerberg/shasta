// Main program for the Shasta static executable.
// The static executable provides
// basic functionality and reduced performance. 
// For full functionality use the shared library built
// under directory src.


#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



int main(int argumentCount, const char** arguments)
{

    Assembler assembler("Data", 2*1024*1024, true);
    return 0;
}
