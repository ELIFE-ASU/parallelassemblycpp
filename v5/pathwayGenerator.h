/**
 * @brief One immutable duplication retained in the best pathway witness.
 *
 * Path steps live only on the active DFS stack and in the winning witness;
 * they are deliberately independent of assembly states and cache entries.
 */
struct assemblyPathStep
{
    EdgeMask match;
    EdgeMask duplicate;
};

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
 * @param step Matching pair to output
 * @param ofs Output file stream
 */
void printMatching(const assemblyPathStep &step, ofstream &ofs)
{
    ofs << "{\"Left\":";
    printMaskAsEdgeList(step.match, ofs);
    ofs << ",\"Right\":";
    printMaskAsEdgeList(step.duplicate, ofs);
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
 * @param pathway Immutable root-to-leaf duplication witness
 * @param removedEdges Edges removed during preprocessing
 * @return true if the complete pathway output was written successfully.
 */
bool recoverPathway2(
    const vector<assemblyPathStep> &pathway,
    const vector<edgeL> &removedEdges
)
{
    EdgeMask allTakenEdges = 0;
    for (const assemblyPathStep &step : pathway)
        allTakenEdges |= step.duplicate;
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
    for (size_t i = 0; i < pathway.size(); i++)
    {
        printMatching(pathway[i], ofs);
        if (i < pathway.size() - 1) ofs << ",\n";
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
