#include <istream>

/**
 * @brief Graph I/O
 * 
 * @param ifs input file
 * @param mg output molGraph
 */
void graphio(std::istream &ifs, molGraph &mg)
{
    string s;
    string graphName;
    int graphSize, a, b;
    vector<pii> edgeList;
    getline(ifs, graphName);
    istringstream nameLine(graphName);
    nameLine >> graphName;
    if (verbose) cout << "Name of graph is: " << graphName << '\n';
    getline(ifs, s);
    istringstream iss1(s);
    iss1 >> graphSize;
    getline(ifs, s);
    istringstream iss(s);
    while (iss >> a >> b)
    {
        edgeList.push_back(pii(a, b));
    }
    getline(ifs, s);
    istringstream iss2(s);
    for (int i = 0; i < graphSize; i++)
    {
        iss2 >> s;
        mg.addAtom(s);
    }
    getline(ifs, s);
    istringstream iss3(s);
    for (const pii &edge : edgeList)
    {
        iss3 >> a;
        mg.addBond(edge.first - 1, edge.second - 1, a);
    }
    if (verbose) mg.printToCout();
}
