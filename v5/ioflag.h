/**
 * @brief runTime flag
 * 
 */
void owTime(string& _runtimeMax)
{
    runTimeMax = atoll(_runtimeMax.c_str());
}

/**
 * @brief subgraph enumeration limit flag
 * 
 */
void owEnumMax(string& _enum_max)
{
    ENUM_MAX = stoi(_enum_max);
}

/**
 * @brief pathway flag
 * 
 */
void owPathway(string& _isPathway)
{
    isPathway = stoi(_isPathway);
}

/**
 * @brief hydrogen removal flag
 * 
 */
void owRemoveHydrogens(string &_removeHydrogens)
{
    removeHydrogens = stoi(_removeHydrogens);
}

/**
 * @brief disjoint graph JAI compensation flag
 * 
 */
void owDisjointCompensate(string &_removeHydrogens)
{
    disjointCompensation = stoi(_removeHydrogens);
}

/**
 * @brief memory test flag
 * 
 */
void owMemTest(string &_memTest)
{
    memTest = stoi(_memTest);
}

/**
 * @brief write intermediate MAs flag
 * 
 */
void owIntermediateMAs(string &_intermediateMA)
{
    writeIntermediateMAs = stoi(_intermediateMA);
}

/// For parsing flags
std::unordered_map<string, void(*)(string&)> fptrTable;

/**
 * @brief For parsing flags
 * 
 */
void fillFptrTable()
{
    void (*f)(string&);
    f = &owTime;
    fptrTable[string("runTime")] = f;
    f = &owEnumMax;
    fptrTable[string("enumMax")] = f;
    f = &owPathway;
    fptrTable[string("pathway")] = f;
    f = &owRemoveHydrogens;
    fptrTable[string("removeHydrogens")] = f;
    f = &owDisjointCompensate;
    fptrTable[string("compensateDisjoint")] = f;
    f = &owMemTest;
    fptrTable[string("testMemory")] = f;
    f = &owIntermediateMAs;
    fptrTable[string("writeIntermediateMAs")] = f;
}

/**
 * @brief For parsing flags
 * 
 */
void flagParser(int argc, char** argv)
{
    fillFptrTable();
    vector<string> args(argc);
    for (int i = 0; i < argc; i++)
    {
        args[i] = argv[i];
        if (args[i][0] == '-')
        {
            string flag = args[i].substr(1, args[i].find("=") - 1);
            string par = args[i].substr(args[i].find("=") + 1, args[i].length() - flag.length() + 2);
            if (fptrTable.count(flag))
            {
                void(*f)(string&) = fptrTable[flag];
                f(par);
            }
        }
    }
}