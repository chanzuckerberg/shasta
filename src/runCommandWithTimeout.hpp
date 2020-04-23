#ifndef SHASTA_RUN_COMMAND_WITH_TIMEOUT_HPP
#define SHASTA_RUN_COMMAND_WITH_TIMEOUT_HPP

#include "string.hpp"

namespace shasta {

    void runCommandWithTimeout(

        // The command to run.
        const string& command,

        // The timeout in seconds
        double timeout,

        // On return, true if the command did not finish
        // and was interrupted because of the timeout.
        bool& timeoutTriggered,

        // On return, true if the command was interrrupted by a signal.
        bool& signalOccurred,

        // Return code from the command.
        // Only valid if timeoutTriggered and signalOccurred are both false.
        // Otherwise, set to -1.
        int& returnCode
        );
}


#endif
