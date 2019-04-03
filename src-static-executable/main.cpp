// Main program for the Shasta static executable.
// The static executable provides
// basic functionality and reduced performance. 
// For full functionality use the shared library built
// under directory src.


#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

#include "iostream.hpp"
#include "stdexcept.hpp"



int main(int argumentCount, const char** arguments)
{
    cout << "The Shasta static executable is not yet functional." << endl;
    try {
        Assembler assembler("Data/", 2*1024*1024, true);
    } catch (exception e) {
        cout << e.what() << endl;
        cout << "Terminated after catching a standard exception." << endl;
        return 1;
    } catch (...) {
        cout << "Terminated after catching a non-standard exception." << endl;
        return 2;
    }
    return 0;
}
