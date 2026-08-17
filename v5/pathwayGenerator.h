/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param mask Bitset representing edges to output
 * @param ofs Output file stream
 */
void printMaskAsEdgeList(EdgeMask mask, ofstream &ofs)
{
    ofs << "[";
    bool first = true;
    for (size_t i = 0; i < univEdgeList.size(); i++)
    {
        if (mask[i])
        {
            if (!first) ofs << ',';
            ofs << "[" << univEdgeList[i].a << "," << univEdgeList[i].b << "]";
            first = false;
        }
    }
    ofs << "]";
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param ofs Output file stream
 */
void printMaskAsEdgeList(ofstream &ofs)
{
    ofs << "[";
    bool first = true;
    for (const edgeL &edge : originalEdgeList)
    {
        if (!first) ofs << ',';
        ofs << "[" << edge.a << "," << edge.b << "]";
        first = false;
    }
    ofs << "]";
}

/**
 * @brief Write the JSON representation used for a bond type.
 */
void printBondColour(short type, ofstream &ofs)
{
    switch (type)
    {
        case 0:
            ofs << "error";
            break;
        case 1:
            ofs << "\"single\"";
            break;
        case 2:
            ofs << "\"double\"";
            break;
        case 3:
            ofs << "\"triple\"";
            break;
        default:
            if (type >= 4) ofs << '"' << type << '"';
            break;
    }
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param maskList List of two bitsets (matching pair)
 * @param ofs Output file stream
 */
void printMatching(vector<EdgeMask> &maskList, ofstream &ofs)
{
    ofs << "{\"Left\":";
    printMaskAsEdgeList(maskList[0], ofs);
    ofs << ",\"Right\":";
    printMaskAsEdgeList(maskList[1], ofs);
    ofs << "}";
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param ofs Output file stream
 */
void printOriginalGraph(ofstream &ofs)
{
    ofs << "\"Vertices\": [";
    for (size_t i = 0; i < originalMolecule.mg.size(); i++)
    {
            ofs << i;
            if (i < originalMolecule.mg.size() - 1) ofs << ',';
    }
    ofs << "],\n";
    ofs << "\"Edges\": ";
    printMaskAsEdgeList(ofs);
    ofs << ",\n";
    ofs << "\"VertexColours\": [";
    for (size_t i = 0; i < originalMolecule.mg.size(); i++)
    {
        ofs << "\"" << originalMolecule.mg[i].type << "\"";
        if (i < originalMolecule.mg.size() - 1) ofs << ',';
    }
    ofs << "],\n";
    ofs << "\"EdgeColours\": [";
    for (size_t i = 0; i < originalEdgeList.size(); i++)
    {
        printBondColour(
            originalMolecule.btypeS(originalEdgeList[i].a, originalEdgeList[i].c),
            ofs
        );
        if (i < originalEdgeList.size() - 1) ofs << ',';
    }
    ofs << "]\n";
}


/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param mask Bitset representing target-graph edges to remove
 * @param ofs Output file stream
 */
void printRemnantGraph(EdgeMask mask, ofstream &ofs)
{
    EdgeMask dual = allEdges ^ mask;
    AtomMask remnantAtoms = 0;
    for (size_t i = 0; i < univEdgeList.size(); i++)
    {
        if (dual[i])
        {
            remnantAtoms.set(univEdgeList[i].a);
            remnantAtoms.set(univEdgeList[i].b);
        }
    }
    ofs << "\"Vertices\": [";
    bool first = true;
    for (size_t i = 0; i < targetMolecule.mg.size(); i++)
    {
        if (remnantAtoms[i])
        {
            if (!first) ofs << ',';
            ofs << i;
            first = false;
        }
    }
    ofs << "],\n";
    ofs << "\"Edges\": ";
    printMaskAsEdgeList(dual, ofs);
    ofs << ",\n";
    ofs << "\"VertexColours\": [";
    first = true;
    for (size_t i = 0; i < targetMolecule.mg.size(); i++)
    {
        if (remnantAtoms[i])
        {
            if (!first) ofs << ',';
            ofs << "\"" << targetMolecule.mg[i].type << "\"";
            first = false;
        }
    }
    ofs << "],\n";
    ofs << "\"EdgeColours\": [";
    first = true;
    for (size_t i = 0; i < univEdgeList.size(); i++)
    {
        if (dual[i])
        {
            if (!first) ofs << ',';
            printBondColour(
                targetMolecule.btypeS(univEdgeList[i].a, univEdgeList[i].c),
                ofs
            );
            first = false;
        }
    }
    ofs << "]\n";
}

/**
 * @brief Pathway reconstruction function, which outputs the original graph, remnants and duplicates to a file
 *
 * Outputs the pathway to the file named by moleculeName, normally INPUTPathway.
 *
 * @param removedEdges Edges removed during preprocessing
 * @return true if the complete pathway output was written successfully.
 */
bool recoverPathway2(vector<edgeL> &removedEdges)
{
    vector<assemblyPath*> minPath;
    for (assemblyPath *current = minAssemblyPath; current != nullptr; current = current->parent)
    {
        minPath.push_back(current);
    }
    reverse(minPath.begin(), minPath.end());
    vector<vector<EdgeMask> > maskList(graphHashMap.size());
    for (auto it = graphHashMap.begin(); it != graphHashMap.end(); ++it)
    {
        maskList[it->second.first].resize(it->second.second + 1);
    }
    for (auto it = bitsetHashTable.begin(); it != bitsetHashTable.end(); ++it)
    {
        maskList[it->second.first][it->second.second] = it->first;
    }
    EdgeMask allTakenEdges = 0;
    vector<vector<EdgeMask>> matchings;
    for (size_t i = 1; i < minPath.size(); i++)
    {
        EdgeMask mask = maskList[minPath[i]->key[0]][minPath[i]->match],
        duplicate = maskList[minPath[i]->key[0]][minPath[i]->duplicate];
        allTakenEdges |= duplicate;
        matchings.push_back({mask, duplicate});
    }
    ofstream ofs(moleculeName);
    if (!ofs.is_open())
    {
        cerr << "error: could not open output file '" << moleculeName << "'\n";
        return false;
    }

    ofs << "{\n";
    ofs << "\"file_graph\":[\n";
    ofs << "{\n";
    printOriginalGraph(ofs);
    ofs << "}\n";
    ofs << "],\n";
    ofs << "\"remnant\":[\n";
    ofs << "{\n";
    printRemnantGraph(allTakenEdges, ofs);
    ofs << "}\n";
    ofs << "],\n";
    ofs << "\"duplicates\":[\n";
    for (size_t i = 0; i < matchings.size(); i++)
    {
        printMatching(matchings[i], ofs);
        if (i < matchings.size() - 1) ofs << ",\n";
    }
    ofs << "\n],\n";
    ofs << "\"removed_edges\":[";    
    for (size_t i = 0; i < removedEdges.size(); i++)
    {
        ofs << "[" << removedEdges[i].a << "," << removedEdges[i].b << "]";
        if (i < removedEdges.size() - 1) ofs << ',';
    }
    ofs << "]\n";
    ofs << "}\n";

    ofs.close();
    if (!ofs)
    {
        cerr << "error: could not write output file '" << moleculeName << "'\n";
        return false;
    }

    return true;
}
