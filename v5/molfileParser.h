/**
 * @brief Takes a molfile and turns it into a molGraph and returns the name of the molecule
 *
 * @param molfile the .mol input
 * @param mg the molGraph output of the function
 * @return string name of molecule
 */
string molfileParser(ifstream &molfile, molGraph &mg)
{
    string name, currLine;

    getline(molfile, name);
    for (int i = 0; i < 2; i++) {getline(molfile, currLine);}
    
    getline(molfile, currLine);
    const int atomCount = stoi(currLine.substr(0, 3));
    const int bondCount = stoi(currLine.substr(3, 3));
    cout << "Detecting " << atomCount << " atoms and " << bondCount << " bonds\n";

    for (int i = 0; i < atomCount; i++)
    {
        getline(molfile, currLine);
        istringstream atomLine(currLine);
        double ignoredCoordinate;
        for (int coordinate = 0; coordinate < 3; coordinate++)
        {
            atomLine >> ignoredCoordinate;
        }
        string atomType;
        atomLine >> atomType;
        mg.addAtom(atomType);
    }
    for (int i = 0; i < bondCount; i++)
    {
        getline(molfile, currLine);
        const int atomA = stoi(currLine.substr(0, 3));
        const int atomB = stoi(currLine.substr(3, 3));
        const int bondOrder = stoi(currLine.substr(6, 3));
        mg.addBond(atomA - 1, atomB - 1, bondOrder);
    }
    if (removeHydrogens)
    {
        for (size_t i = 0; i < mg.mg.size(); i++)
        {
            if (mg.atype(i) == "H") mg.removeAtom(i);
        }
        mg.removeAndCollapse();
    }
    mg.printToCout();
    return name;
}
