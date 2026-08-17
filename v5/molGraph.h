/**
 * @brief Bond struct for molGraph
 */
struct bond
{
    short n;
    short type;
    bond(short _n, short _type): n(_n), type(_type){}
};

/**
 * @brief Atom struct for molGraph
 */
struct atom
{
    string type;
    vector<bond> list;

    atom(string _type): type(_type){}

};

/**
 * @brief Primary graph data structure used in assemblyCpp
 */
struct molGraph
{
    /**
     * @brief Vector of atoms representing nodes in the graph.
     */
    vector<atom> mg;
    /**
     * @brief Total number of bonds (edges) in the graph.
     */
    int totalBonds = 0;

    /**
     * @brief Use this function to add atoms/nodes
     * @param _type Type is atom type/node labelling.
     */
    void addAtom(string &_type)
    {
        atom a(_type);
        mg.push_back(a);
    }

    /**
     * @brief Use this function to add bonds/edges.
     * @param a Index of first atom/node
     * @param b Index of second atom/node
     * @param type Type is bond order/edge labelling.
     */
    void addBond(int a, int b, short type)
    {
        bond b1(b, type), b2(a, type);
        mg[a].list.push_back(b1);
        mg[b].list.push_back(b2);
        totalBonds++;
    }
    
    /**
     * @brief Print the graph information to cout
     */
    void printToCout()
    {
        cout << "There are " << mg.size() << " atoms in the molecule-graph\n";
        for (size_t i = 0; i < mg.size(); i++)
        {
            cout << "Atom " << i + 1 << " is of type " << mg[i].type << " and adjacent to atoms ";
            for (size_t j = 0; j < degree(i); j++)
            {
                cout << elem(i, j) + 1 << " with bond order " << btypeS(i, j) << ", ";
            }
            cout << '\n';
        }
    }

    /**
     * @brief Get the degree (number of bonds) of atom at index x
     */
    size_t degree(int x)
    {
        return mg[x].list.size();
    }

    short elem(size_t a, size_t b)
    {
        return mg[a].list[b].n;
    }

    /**
     * @brief Get atom type for index i
     */
    string atype(size_t i) {return mg[i].type;}

    /**
     * @brief Get bond type as char
     * @param a Index of first atom/node
     * @param b Index of bond
     */
    char btype(size_t a, size_t b)
    {
        return static_cast<char>(mg[a].list[b].type);
    }

    /**
     * @brief Get bond type as short
     * @param a Index of first atom/node
     * @param b Index of bond
     */
    short btypeS(size_t a, size_t b)
    {
        return mg[a].list[b].type;
    }

private:
    /**
     * @brief Rebuild the graph while dropping removed atoms and zero-order bonds.
     */
    void rebuild(bool removeMarkedAtoms)
    {
        const size_t removed = mg.size();
        vector<size_t> reverseMap(mg.size(), removed);
        molGraph output;

        for (size_t i = 0; i < mg.size(); i++)
        {
            if (removeMarkedAtoms && mg[i].type == "COLLAPSE") continue;
            reverseMap[i] = output.mg.size();
            output.addAtom(mg[i].type);
        }

        for (size_t source = 0; source < mg.size(); source++)
        {
            if (reverseMap[source] == removed) continue;

            for (size_t bondIndex = 0; bondIndex < degree(source); bondIndex++)
            {
                const size_t target = static_cast<size_t>(elem(source, bondIndex));
                if (
                    reverseMap[target] != removed &&
                    btypeS(source, bondIndex) != 0 &&
                    reverseMap[source] < reverseMap[target]
                )
                {
                    output.addBond(
                        static_cast<int>(reverseMap[source]),
                        static_cast<int>(reverseMap[target]),
                        btype(source, bondIndex)
                    );
                }
            }
        }

        *this = std::move(output);
    }

public:

    /**
     * @brief Remove all bonds of order 0. Used in preprocessing.
     */
    void collapse()
    {
        rebuild(false);
    }

    /**
     * @brief For explicit hydrogen removal
     * 
     */
    void removeAtom(size_t i)
    {
        if (i >= mg.size()) return;
        mg[i].type = "COLLAPSE";
    }

    /**
     * @brief For explicit hydrogen removal
     * 
     */
    void removeAndCollapse()
    {
        rebuild(true);
    }

    /**
     * @brief Turns molGraph (adjacency list) into equivalent edgelist
     * @return std::vector<edgeL>
     */
    vector<edgeL> writeEdgeList()
    {
        vector<edgeL> out;
        for (size_t i = 0; i < mg.size(); i++)
        {
            for (size_t j = 0; j < degree(i); j++)
            {
                short k = elem(i, j);
                if (i < static_cast<size_t>(k))
                {
                    short source = static_cast<short>(i);
                    short bondIndex = static_cast<short>(j);
                    edgeL t(source, k, bondIndex);
                    out.push_back(t);
                }
            }
        }
        return out;
    }
    
    /**
     * @brief For preprocessing, writes edgeList as hash map to detect duplicated bonds
     */
    void writeEdgeList(std::unordered_map<string, pair<int, edgeL> > &ht)
    {
        for (size_t i = 0; i < mg.size(); i++)
        {
            for (size_t j = 0; j < degree(i); j++)
            {
                short k = elem(i, j);
                if (i < static_cast<size_t>(k))
                {
                    string is = atype(i), ks = atype(k), out;
                    if (is < ks) out = is + btype(i, j) + ks;
                    else out = ks + btype(i, j) + is;
                    short source = static_cast<short>(i);
                    short bondIndex = static_cast<short>(j);
                    edgeL t(source, k, bondIndex);
                    auto [entry, inserted] = ht.try_emplace(out, 1, t);
                    if (!inserted) entry->second.first++;
                }
            }
        }
    }

    /**
     * @brief Used in preprocessing, removes edges in edgelist
     */
    void negativeEdgeCollapse(vector<edgeL> &edgeList)
    {
        for (size_t i = 0; i < edgeList.size(); i++)
        {
            edgeL &el = edgeList[i];
            mg[el.a].list[el.c].type = 0;
        }
        collapse();
    }

    /**
     * @brief For compensating for disjoint fragments in the JAI
     * 
     */
    void disjointFragmentsR(vb & visited, int n)
    {
        if (visited[n]) return;
        visited[n] = 1;
        for (size_t i = 0; i < mg[n].list.size(); i++) disjointFragmentsR(visited, elem(n, i));
    }

    /**
     * @brief For compensating for disjoint fragments in the JAI
     * 
     * @return int number of disjoint fragments
     */
    int disjointFragments()
    {
        vb visited(mg.size(), 0); int count = 0;
        for (size_t i = 0; i < mg.size(); i++)
        {
            if (!visited[i])
            {
                count++;
                disjointFragmentsR(visited, i);
            }
        }
        return count;
    }
};

/**
 * @brief construct new molGraph from input molGraph and boolean edgelist
 *
 * @param mg Input molgraph
 * @param edgeList Corresponding edge list
 * @param mask Input boolean edgelist
 * @param isCyclic Is the graph cyclic, needed for hashing
 * @return molGraph
 */
molGraph constructFromEdgeList(molGraph &mg, vector<edgeL> &edgeList, 
    standardBitset &mask, bool &isCyclic)
    {
        disjointSet u(mg.mg.size());
        molGraph output;
        std::unordered_map<int, int> ht;
        isCyclic = 0;
        for (size_t i = 0; i < edgeList.size(); i++)
        {
            if (mask[i] != 0)
            {
                int a = edgeList[i].a, b = edgeList[i].b, c = 0;
                if (ht.count(a) == 0 && ht.count(b) == 0)
                {
                    size_t x = ht.size();
                    ht[a] = x;
                    output.addAtom(mg.mg[a].type);
                    size_t y = ht.size();
                    ht[b] = y;
                    output.addAtom(mg.mg[b].type);
                    u.insert(x, x);
                    u.insert(y, x);
                }
                else
                {
                    if (ht.count(a) == 0)
                    {
                        size_t x = ht.size();
                        ht[a] = x;
                        output.addAtom(mg.mg[a].type);
                        u.insert(x, ht[b]);
                    }
                    else c++;
                    if (ht.count(b) == 0)
                    {
                        size_t x = ht.size();
                        ht[b] = x;
                        output.addAtom(mg.mg[b].type);
                        u.insert(x, ht[a]);
                    }
                    else c++;
                }
                int a2 = ht[a], b2 = ht[b];
                output.addBond(a2, b2, mg.btype(a, edgeList[i].c));
                if (c == 2)
                {
                    isCyclic |= u.merge(a2, b2);
                }
            }
        }
        return output;
    }

/**
 * @brief Preprocesses the graph by removing all unique edges for the pathway algorithm
 * @param mg The input molGraph
 * @param writeback Edges removed during preprocessing
 * @return molGraph (the final output)
 */
molGraph preprocessWriteback(molGraph &mg, vector<edgeL> &writeback)
{
    std::unordered_map<string, pair<int, edgeL> > ht;
    molGraph out = mg;
    mg.writeEdgeList(ht);
    vector<edgeL> v;
    for (auto it = ht.begin(); it != ht.end(); ++it)
    {
        if (it->second.first == 1)
        {
            v.push_back(it->second.second);
        }
    }
    out.negativeEdgeCollapse(v);
    writeback = v;
    return out;
}

/// Global variable for the molGraph before and after preprocessing
molGraph originalMolecule, targetMolecule;

struct ufdsMaskWorkspace
{
    ufdsSplit sets;
    vi uniques;
    vector<standardBitset> components;

    ufdsMaskWorkspace(size_t atomCount, size_t edgeCount)
    {
        sets.elements.reserve(atomCount);
        sets.extraVals.reserve(edgeCount);
        uniques.reserve(atomCount);
        components.reserve(atomCount);
    }
};

/**
 * @brief Calls the disjoint-set data structure for the fragmentation function. See Seet et al. section 4.5 for details
 *
 * @param mask Target bitset as input
 * @param maskList List of disjoint bitsets returned
 * @param workspace Reusable disjoint-set and component buffers. Its component
 * buffer must not alias maskList.
 */
void ufdsMaskConstructWithWorkspace(
    standardBitset &mask,
    vector<standardBitset> &maskList,
    ufdsMaskWorkspace &workspace
)
{
    vector<edgeL> &edgeList = univEdgeList;
    ufdsSplit &u = workspace.sets;
    u.reset(targetMolecule.mg.size());
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        if (mask[i] != 0)
        {
            int a = edgeList[i].a, b = edgeList[i].b;
            if (u.elements[a].parent == -1 && u.elements[b].parent == -1)
            {
                u.doubleInsert(b, a, i);
            }
            else
            {
                if (u.elements[a].parent == -1)
                {
                    u.insert(a, b, i);
                }
                else if (u.elements[b].parent == -1)
                {
                    u.insert(b, a, i);
                }
                else
                {
                    u.merge(a, b, i);
                }
            }
        }
    }
    u.splitWithBuffers(maskList, workspace.uniques, workspace.components);
}
