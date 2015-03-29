#include "SignalHandlers.h"
#include <signal.h>
#include <iostream>

using namespace std;
void (*_cleanup)() = NULL;

void HandleSignals(int s){
    switch (s)
    {
        case SIGHUP:
        case SIGTERM:
        case SIGSEGV:
        case SIGABRT:
        case SIGFPE:
        case SIGILL:
        case SIGINT:
            try
            {
            	cerr << "Handling signal: " << s << endl << flush;
                if (_cleanup != NULL)
                    (*_cleanup)();
            	cerr << "done!" << endl << flush;
            }
            catch (...)
            {}
            break;
        default:
            break;
    }
}

void ConfigureSignalHandlers(void (*cleanup)())
{
    _cleanup = cleanup;

    struct sigaction sigHandler;
    sigHandler.sa_handler = HandleSignals;
    sigemptyset(&sigHandler.sa_mask);
    sigHandler.sa_flags = 0;

    sigaction(SIGHUP, &sigHandler, NULL);
    sigaction(SIGTERM, &sigHandler, NULL);
    sigaction(SIGSEGV, &sigHandler, NULL);
    sigaction(SIGABRT, &sigHandler, NULL);
    sigaction(SIGFPE, &sigHandler, NULL);
    sigaction(SIGILL, &sigHandler, NULL);
    sigaction(SIGINT, &sigHandler, NULL);
};

void RemoveSignalHandlers()
{
    _cleanup = NULL;    

    struct sigaction sigHandler;
    sigHandler.sa_handler = SIG_DFL;
    sigemptyset(&sigHandler.sa_mask);
    sigHandler.sa_flags = 0;

    sigaction(SIGHUP, &sigHandler, NULL);
    sigaction(SIGTERM, &sigHandler, NULL);
    sigaction(SIGSEGV, &sigHandler, NULL);
    sigaction(SIGABRT, &sigHandler, NULL);
    sigaction(SIGFPE, &sigHandler, NULL);
    sigaction(SIGILL, &sigHandler, NULL);
    sigaction(SIGINT, &sigHandler, NULL);
};
