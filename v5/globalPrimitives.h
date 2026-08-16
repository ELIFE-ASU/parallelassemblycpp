template <typename T1, typename T2, typename T3>
struct triple
{
    T1 a; T2 b; T3 c;
    triple(){}
    triple(T1 &_a, T2 &_b, T3 &_c): a(_a), b(_b), c(_c) {}
};

int DISCHARGE_FREQUENCY = 30000000, ENUM_MAX = 50000000;
int minAIfound = -1;
int recursiveCount = 1;
int assemblyIx = 0;

std::unordered_map<string, int> atypeHash;
vector<double> coords;
string moleculeName;
standardBitset allEdges;
volatile bool interruptFlag = false;
clock_t startTime = 0;
unsigned long long runTimeMax = 18446744073709551615;

typedef triple<int, int, int> iii;
typedef triple<short, short, short> edgeL;
unsigned int totalBonds = 0;
vector<edgeL> removedEdges;
vector<edgeL> originalEdgeList, univEdgeList;

/// Hash table for edgelists for pathway algorithm
std::unordered_map<standardBitset, pii> bitsetHashTable;

bool isPathway = 1, removeHydrogens = 1, disjointCompensation = 0, memTest = 0, writeIntermediateMAs = 0;
int disjointFragments = 1;
vector<pair<unsigned long long, int> > intermediateMAs;