#include "runCommandWithTimeout.hpp"
#include "platformDependent.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

#include <cstdlib>



void shasta::runCommandWithTimeout(

    // The command to run.
    const string& command,

    // The timeout in seconds
    double timeout,

    // On return, true if the command did not finish
    // and was interrupted because of the timeout.
    bool& timeoutTriggered,

    // On return, true if the command was interrupted by a signal.
    bool& signalOccurred,

    // Return code from the command.
    // Only valid if timeoutTriggered and signalOccurred are both false.
    // Otherwise, set to -1.
    int& returnCode
    )
{
    timeoutTriggered = false;
    signalOccurred = false;
    returnCode = -1;

    // Handle a negative or zero timeout.
    // You can't get much done in negative or zero time.
    if(timeout <= 0.) {
        timeoutTriggered = true;
        return;
    }

    // Run the command using the timeout command.
    const string commandWithTimeout =
        timeoutCommand() + " " + to_string(timeout) + " " + command;
    const int commandStatus = std::system(commandWithTimeout.c_str());


    if(WIFEXITED(commandStatus)) {

        // The timeout command completed without being interrupted by a signal.

        // Extract the return code.
        returnCode = WEXITSTATUS(commandStatus);

        if(returnCode == 124) {

            // Interrupted by timeout.
            timeoutTriggered = true;
            returnCode = -1;
            return;

        } else {

            // Normal termination.
            return;
        }

    } else if(WIFSIGNALED(commandStatus)) {

        // Interrupted by a signal.
        // const int signalNumber = WTERMSIG(commandStatus);
        signalOccurred = true;
        return;

    } else {
        SHASTA_ASSERT(0);
    }
}



