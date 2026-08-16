#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
using namespace boost;


typedef property<edge_name_t, char> bond_vf2;
typedef property<vertex_name_t, string, property<vertex_index_t, int> > atom_vf2;

typedef adjacency_list <vecS, vecS, undirectedS, atom_vf2, bond_vf2> molGraphBoost;
typedef property_map<molGraphBoost, vertex_name_t >::type vertex_name_map_t;
typedef property_map_equivalent< vertex_name_map_t, vertex_name_map_t >
    vertex_comp_t;
typedef property_map<molGraphBoost, edge_name_t>::type edge_name_map_t;
typedef property_map_equivalent<edge_name_map_t, edge_name_map_t> edge_comp_t;

/**
 * @brief Converts a boolean edgelist based on a standard molGraph and edgeList into a Boost molGraph
 *
 * @param mg The standard molGraph
 * @param edgeList The standard edgeList
 * @param mask The input boolean edgelist
 * @return molGraphBoost the output
 */
molGraphBoost edgelistToBoost(molGraph &mg, vector<edgeL> &edgeList, const standardBitset &mask)
{
    molGraphBoost output;
    std::unordered_map<int, int> ht;
    for (size_t i = 0; i < edgeList.size(); i++)
    {
        if (mask[i])
        {
            int a = edgeList[i].a, b = edgeList[i].b;
            if (ht.count(a) == 0)
            {
                size_t x = ht.size();
                ht[a] = x;
                string c = mg.mg[a].type;
                add_vertex(atom_vf2(c), output);
            }
            if (ht.count(b) == 0)
            {
                size_t x = ht.size();
                ht[b] = x;
                string c = mg.mg[b].type;
                add_vertex(atom_vf2(c), output);
            }
            int a2 = ht[a], b2 = ht[b];
            add_edge(a2, b2, mg.btype(a, edgeList[i].c), output);
        }
    }
    return output;
}

/**
 * @brief Halt as soon as any isomorphism is found
 *
 */
struct halting_callback {
    template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
    bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) const {
        return false;
    }
};

/**
 * @brief vf2 graph isomorphism caller
 *
 * @param mmg First graph to be compared
 * @param tmg Second graph to be compared
 * @return true if the two graphs are isomorphic
 * @return false otherwise
 */
bool vf2GraphIso(molGraphBoost &mmg, molGraphBoost &tmg)
{
    vertex_comp_t vc =
    make_property_map_equivalent(get(vertex_name, mmg), get(vertex_name, tmg));
    edge_comp_t ec =
    make_property_map_equivalent(get(edge_name, mmg), get(edge_name, tmg));
    halting_callback callback;
    return vf2_graph_iso(mmg, tmg, callback, vertex_order_by_mult(mmg), edges_equivalent(ec).vertices_equivalent(vc));
}