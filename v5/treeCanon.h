/**
 * @brief Subroutine of tree canonisation function
 *
 * @param mg Molecular graph (assumed to be a tree)
 * @param curr Current node index
 * @param parent Parent node index
 * @param type Bond type leading to this node
 * @return Canonical string for subtree rooted at curr
 */
string treeCanonRecursive(molGraph &mg, int curr, int parent, char type)
{
    string atomName(mg.mg[curr].type), bondName;
    if (atomName == "X") return "";
    if (parent != -1) bondName = type;
    else bondName = "$";
    string name = "(" + bondName + atomName;
    
    if (mg.degree(curr) == 0) {name += ')'; 
    return name;}
    vector<string> stringArray;
    vector<bond> &v = mg.mg[curr].list;
    for (size_t i = 0; i < v.size(); i++)
    {
        if (v[i].n != parent) stringArray.push_back(treeCanonRecursive(mg, v[i].n, curr, v[i].type));
    }
    sort(stringArray.begin(), stringArray.end());
    for (unsigned int i = 0; i < stringArray.size(); i++) name += stringArray[i];
    name += ')';
    return name;
}

/**
 * @brief Subroutine of centroid finding function
 *
 * @param mg Molecular graph
 * @param weight Vector of subtree weights
 * @param curr Current node index
 * @param prev Previous node index
 * @return Subtree size rooted at curr
 */
int centroidDFS(molGraph &mg, vector<int> &weight, int curr, int prev)
{
    if (mg.mg[curr].type == "X") return 0;
    vector<bond> &v = mg.mg[curr].list;
    for (size_t i = 0; i < v.size(); i++)
    {  
        if (v[i].n != prev) weight[curr] += centroidDFS(mg, weight, v[i].n, curr);
    }
    return weight[curr];
}

/**
 * @brief Subroutine of tree canonisation function, finds centroid
 *
 * @param mg Molecular graph
 * @param currRoot Starting index to search for centroid
 * @return Pair of centroid(s) (second = -1 if only one)
 */
pii centroid(molGraph &mg, int currRoot)
{
    float size = mg.mg.size() + 0.0;
    vector<bool> visited((int)size, 0);
    vector<int> weight((int)size, 1);
    bool isCentroid = 0; int prevRoot = -1;
    centroidDFS(mg, weight, currRoot, -1);
    bool has2Centroids = 0;
    while (!isCentroid)
    {
        vector<bond> &v = mg.mg[currRoot].list;
        bool overweight = 0;
        short j = currRoot;
        for (size_t i = 0; i < v.size(); i++)
        {
            j = v[i].n;
            if (j != prevRoot && weight[j] >= size/2)
            {
                if (weight[j] == size/2) {has2Centroids = 1;}
                overweight = 1; break;
            }
        }
        if (!overweight) isCentroid = 1;
        else {prevRoot = currRoot; currRoot = j;}
    }
    pii p(currRoot, -1);
    if (has2Centroids) p.second = prevRoot;
    return p;
}

/**
 * @brief Function to find a canonical string for a tree, using standard tree canonisation algorithm
 *
 * @param mg Acyclic molGraph
 * @param n Starting index
 * @return string (canonical form)
 */
string centroidTreeCanon(molGraph &mg, int n)
{
    if (mg.mg.size() == 0) return "";
    string canonicalForm;
    pii centroids = centroid(mg, n);
    vector<string> stringArray(2, "");
    stringArray[0] = treeCanonRecursive(mg, centroids.first, -1, ' ');
    if (centroids.second != -1)
    {
        stringArray[1] = treeCanonRecursive(mg, centroids.second, -1, ' ');
        sort(stringArray.begin(), stringArray.end());
        canonicalForm = stringArray[0] + stringArray[1];
    }
    else canonicalForm = stringArray[0];
    return canonicalForm;
}

