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
        // Signal handlers may only perform async-signal-safe work. The search
        // observes this flag and unwinds normally so enabled outputs are
        // written and flushed by the regular control flow.
        userInterruptReceived = 1;
        interruptFlag = 1;
    }

    void disableInterruptHandler() {
        signal(SIGINT, SIG_IGN);
    }
#endif
