#pragma once

#ifdef _WIN32
    BOOL WINAPI CtrlHandler(DWORD fdwCtrlType) {
        switch (fdwCtrlType) {
            case CTRL_C_EVENT:
                userInterruptReceived.store(true);
                interruptFlag.store(true);
                return TRUE;


            default:
                return FALSE;
        }
    }

    void disableInterruptHandler() {
        // A null handler with TRUE makes the process ignore new Ctrl-C events
        // while final output streams are flushed.
        SetConsoleCtrlHandler(nullptr, TRUE);
    }
#else
    void signalHandler(int) {
        // These atomics are guaranteed lock-free and signal-safe. The search
        // observes the flags and unwinds normally so enabled outputs are
        // written and flushed by the regular control flow.
        userInterruptReceived.store(true);
        interruptFlag.store(true);
    }

    void disableInterruptHandler() {
        signal(SIGINT, SIG_IGN);
    }
#endif
