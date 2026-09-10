#pragma once

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
 * @param outputStream Output file stream
 */
void printMaskAsEdgeList(EdgeMask mask, ofstream &outputStream)
{
    const vector<MoleculeEdge> &edgeList = searchUniverseEdgeList();
    outputStream << "[";
    bool first = true;
    for (size_t edgeIndex = 0; edgeIndex < edgeList.size(); edgeIndex++)
    {
        if (mask[edgeIndex])
        {
            if (!first) outputStream << ',';
            outputStream
                << "[" << edgeList[edgeIndex].sourceAtomIndex << ","
                << edgeList[edgeIndex].targetAtomIndex << "]";
            first = false;
        }
    }
    outputStream << "]";
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param outputStream Output file stream
 */
void printMaskAsEdgeList(ofstream &outputStream)
{
    outputStream << "[";
    bool first = true;
    for (const MoleculeEdge &edge : originalEdgeList)
    {
        if (!first) outputStream << ',';
        outputStream << "[" << edge.sourceAtomIndex << ","
            << edge.targetAtomIndex << "]";
        first = false;
    }
    outputStream << "]";
}

/**
 * @brief Write one JSON string, escaping every required ASCII character.
 */
void printJsonString(const string &value, ostream &output)
{
    static constexpr char hexDigits[] = "0123456789ABCDEF";
    output.put('"');
    for (const unsigned char character : value)
    {
        switch (character)
        {
            case '"': output << "\\\""; break;
            case '\\': output << "\\\\"; break;
            case '\b': output << "\\b"; break;
            case '\f': output << "\\f"; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default:
                if (character < 0x20)
                {
                    output << "\\u00"
                           << hexDigits[character >> 4]
                           << hexDigits[character & 0x0f];
                }
                else output.put(static_cast<char>(character));
                break;
        }
    }
    output.put('"');
}

/**
 * @brief Write the JSON representation used for a bond type.
 */
void printBondColour(short type, ostream &output)
{
    switch (type)
    {
        case 1:
            output << "\"single\"";
            break;
        case 2:
            output << "\"double\"";
            break;
        case 3:
            output << "\"triple\"";
            break;
        default:
            if (type >= 4) output << '"' << type << '"';
            else output << "\"error\"";
            break;
    }
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param step Matching pair to output
 * @param outputStream Output file stream
 */
void printMatching(const assemblyPathStep &step, ofstream &outputStream)
{
    outputStream << "{\"Left\":";
    printMaskAsEdgeList(step.match, outputStream);
    outputStream << ",\"Right\":";
    printMaskAsEdgeList(step.duplicate, outputStream);
    outputStream << "}";
}

/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param outputStream Output file stream
 */
void printOriginalGraph(ofstream &outputStream)
{
    outputStream << "\"Vertices\": [";
    for (
        size_t atomIndex = 0;
        atomIndex < originalMolecule.atoms.size();
        atomIndex++
    )
    {
            outputStream << atomIndex;
            if (atomIndex < originalMolecule.atoms.size() - 1)
                outputStream << ',';
    }
    outputStream << "],\n";
    outputStream << "\"Edges\": ";
    printMaskAsEdgeList(outputStream);
    outputStream << ",\n";
    outputStream << "\"VertexColours\": [";
    for (
        size_t atomIndex = 0;
        atomIndex < originalMolecule.atoms.size();
        atomIndex++
    )
    {
        printJsonString(
            originalMolecule.atoms[atomIndex].atomType,
            outputStream
        );
        if (atomIndex < originalMolecule.atoms.size() - 1)
            outputStream << ',';
    }
    outputStream << "],\n";
    outputStream << "\"EdgeColours\": [";
    for (
        size_t edgeIndex = 0;
        edgeIndex < originalEdgeList.size();
        edgeIndex++
    )
    {
        printBondColour(
            originalMolecule.bondType(
                originalEdgeList[edgeIndex].sourceAtomIndex,
                originalEdgeList[edgeIndex].sourceBondIndex
            ),
            outputStream
        );
        if (edgeIndex < originalEdgeList.size() - 1) outputStream << ',';
    }
    outputStream << "]\n";
}


/**
 * @brief Subroutine of the pathway reconstruction function
 *
 * @param mask Bitset representing target-graph edges to remove
 * @param outputStream Output file stream
 */
void printRemnantGraph(EdgeMask mask, ofstream &outputStream)
{
    const vector<MoleculeEdge> &edgeList = searchUniverseEdgeList();
    const molGraph &molecule = searchTargetMolecule();
    EdgeMask dual = allEdges ^ mask;
    AtomMask remnantAtoms = 0;
    for (size_t edgeIndex = 0; edgeIndex < edgeList.size(); edgeIndex++)
    {
        if (dual[edgeIndex])
        {
            remnantAtoms.set(edgeList[edgeIndex].sourceAtomIndex);
            remnantAtoms.set(edgeList[edgeIndex].targetAtomIndex);
        }
    }
    outputStream << "\"Vertices\": [";
    bool first = true;
    for (size_t atomIndex = 0; atomIndex < molecule.atoms.size(); atomIndex++)
    {
        if (remnantAtoms[atomIndex])
        {
            if (!first) outputStream << ',';
            outputStream << atomIndex;
            first = false;
        }
    }
    outputStream << "],\n";
    outputStream << "\"Edges\": ";
    printMaskAsEdgeList(dual, outputStream);
    outputStream << ",\n";
    outputStream << "\"VertexColours\": [";
    first = true;
    for (size_t atomIndex = 0; atomIndex < molecule.atoms.size(); atomIndex++)
    {
        if (remnantAtoms[atomIndex])
        {
            if (!first) outputStream << ',';
            printJsonString(molecule.atoms[atomIndex].atomType, outputStream);
            first = false;
        }
    }
    outputStream << "],\n";
    outputStream << "\"EdgeColours\": [";
    first = true;
    for (size_t edgeIndex = 0; edgeIndex < edgeList.size(); edgeIndex++)
    {
        if (dual[edgeIndex])
        {
            if (!first) outputStream << ',';
            printBondColour(
                molecule.bondType(
                    edgeList[edgeIndex].sourceAtomIndex,
                    edgeList[edgeIndex].sourceBondIndex
                ),
                outputStream
            );
            first = false;
        }
    }
    outputStream << "]\n";
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
    const vector<MoleculeEdge> &removedEdges
)
{
    EdgeMask allTakenEdges = 0;
    for (const assemblyPathStep &step : pathway)
        allTakenEdges |= step.duplicate;
    ofstream outputStream(moleculeName);
    if (!outputStream.is_open())
    {
        cerr << "error: could not open output file '" << moleculeName << "'\n";
        return false;
    }

    outputStream << "{\n";
    outputStream << "\"file_graph\":[\n";
    outputStream << "{\n";
    printOriginalGraph(outputStream);
    outputStream << "}\n";
    outputStream << "],\n";
    outputStream << "\"remnant\":[\n";
    outputStream << "{\n";
    printRemnantGraph(allTakenEdges, outputStream);
    outputStream << "}\n";
    outputStream << "],\n";
    outputStream << "\"duplicates\":[\n";
    for (size_t stepIndex = 0; stepIndex < pathway.size(); stepIndex++)
    {
        printMatching(pathway[stepIndex], outputStream);
        if (stepIndex < pathway.size() - 1) outputStream << ",\n";
    }
    outputStream << "\n],\n";
    outputStream << "\"removed_edges\":[";
    for (
        size_t removedEdgeIndex = 0;
        removedEdgeIndex < removedEdges.size();
        removedEdgeIndex++
    )
    {
        outputStream
            << "[" << removedEdges[removedEdgeIndex].sourceAtomIndex << ","
            << removedEdges[removedEdgeIndex].targetAtomIndex << "]";
        if (removedEdgeIndex < removedEdges.size() - 1) outputStream << ',';
    }
    outputStream << "]\n";
    outputStream << "}\n";

    outputStream.close();
    if (!outputStream)
    {
        cerr << "error: could not write output file '" << moleculeName << "'\n";
        return false;
    }

    return true;
}
