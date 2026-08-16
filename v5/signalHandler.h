#ifdef _WIN32
    BOOL CtrlHandler(DWORD fdwCtrlType) {
        switch (fdwCtrlType) {
            case CTRL_C_EVENT:
                cout << "min AI found so far: " << minAIfound << '\n';
                interruptFlag = true;
                return TRUE;


            default:
                return FALSE;
        }
    }
#else
    void signalHandler(int signum) {
        cout << "Interrupt signal received. Lowest AI found: " <<  minAIfound << '\n';
        recoverPathway2(removedEdges);
        exit(signum);
    }
#endif