#include <indigox/algorithm/access.hpp>
#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

#include <boost/graph/vf2_sub_graph_iso.hpp>

#include <EASTL/vector_map.h>

namespace indigox::algorithm {
  using namespace indigox::graph;

  bool CMGPrintCallback::operator()(const CorrespondenceMap &cmap) {
    std::cout << "Mapping instance " << ++count << ":\n";
    for (auto &ab : cmap) std::cout << ab.first << " -- " << ab.second << "\n";
    std::cout << "\n";
    return true;
  }

  bool CMGPrintCallback::operator()(const CMGVertex &vs, const CMGVertex &vl) {
    return vs.GetIsomorphismMask() == vl.GetIsomorphismMask();
  }

  bool MGPrintCallback::operator()(const CorrespondenceMap &cmap) {
    std::cout << "Mapping instance " << ++count << ":\n";
    for (auto &ab : cmap) std::cout << ab.first << " -- " << ab.second << "\n";
    std::cout << "\n";
    return true;
  }

  template <class V, class E, class S, class D, class VP, class EP, class VIM>
  struct InternalCallback {
    using GraphType = BaseGraph<V, E, S, D, VP, EP>;
    using BoostGraph = typename GraphType::graph_type;
    using Vertex = typename BoostGraph::vertex_descriptor;
    using Edge = typename BoostGraph::edge_descriptor;
    using UserCallback = MappingCallback<V, E, S, D, VP, EP>;
    using VertDescMap = typename GraphType::VertMap::RightType;
    using EdgeDescMap = typename GraphType::EdgeMap::RightType;

    GraphType &small;
    GraphType &large;
    VIM &v_id_map;
    UserCallback &user_callback;
    VertDescMap &vdm_small;
    VertDescMap &vdm_large;
    EdgeDescMap &edm_small;
    EdgeDescMap &edm_large;

    InternalCallback(GraphType &s, GraphType &l, VIM &idxmap, UserCallback &cb)
        : small(s), large(l), v_id_map(idxmap), user_callback(cb),
          vdm_small(access::GetVertexMap(s).right),
          vdm_large(access::GetVertexMap(l).right),
          edm_small(access::GetEdgeMap(s).right),
          edm_large(access::GetEdgeMap(l).right) {}

    template <class CMap1to2, class CMap2to1>
    bool operator()(CMap1to2 one, CMap2to1) {
      eastl::vector_map<V, V> map;
      typename BoostGraph::vertex_iterator begin, end;
      std::tie(begin, end) = boost::vertices(access::GetGraph(small));
      for (; begin != end; ++begin) {
        V vs = vdm_small.at(*begin);
        V vl = vdm_large.at(boost::get(one, *begin));
        map.emplace(vs, vl);
      }
      return user_callback(map);
    }

    virtual bool operator()(Vertex vs, Vertex vl) {
      return user_callback(vdm_small.at(vs), vdm_large.at(vl));
    }
    virtual bool operator()(Edge es, Edge el) {
      return user_callback(edm_small.at(es), edm_large.at(el));
    }
  };

  template <class V, class E, class S, class D, class VP, class EP>
  void
  SubgraphIsomorphismsRunner(graph::BaseGraph<V, E, S, D, VP, EP> &G1,
                             graph::BaseGraph<V, E, S, D, VP, EP> &G2,
                             MappingCallback<V, E, S, D, VP, EP> &callback) {

    using GraphType = graph::BaseGraph<V, E, S, D, VP, EP>;
    using BoostGraph = typename GraphType::graph_type;
    using Vertex = typename GraphType::VertType;
    using VertIter = typename GraphType::VertIter;
    using VertexMap = eastl::vector_map<Vertex, size_t>;

    if (G2.NumVertices() < G1.NumVertices()) {
      SubgraphIsomorphismsRunner(G2, G1, callback);
      return;
    }

    BoostGraph &small = access::GetGraph(G1);
    BoostGraph &large = access::GetGraph(G2);

    // Make the vertex map
    VertexMap vmap;
    size_t i = 0;
    VertIter begin, end;
    for (std::tie(begin, end) = boost::vertices(small); begin != end; ++begin)
      vmap.emplace(*begin, i++);
    i = 0;
    for (std::tie(begin, end) = boost::vertices(large); begin != end; ++begin)
      vmap.emplace(*begin, i++);

    // Determine the order of vertices
    std::tie(begin, end) = boost::vertices(small);
    std::vector<Vertex> v_order(begin, end);
    std::sort(v_order.begin(), v_order.end(), [&small](Vertex v1, Vertex v2) {
      return boost::degree(v1, small) < boost::degree(v2, small);
    });

    // Make the internal callback
    InternalCallback icallback(G1, G2, vmap, callback);
    boost::associative_property_map<VertexMap> propmap(vmap);

    // Run the isomorphism checking
    boost::vf2_subgraph_iso(small, large, icallback, propmap, propmap, v_order,
                            icallback, icallback);
  }

  void SubgraphIsomorphisms(CondensedMolecularGraph &G1,
                            CondensedMolecularGraph &G2, CMGCallback &CB) {
    SubgraphIsomorphismsRunner(G1, G2, CB);
  }

  void SubgraphIsomorphisms(MolecularGraph &G1, MolecularGraph &G2,
                            MGCallback &CB) {
    SubgraphIsomorphismsRunner(G1, G2, CB);
  }
} // namespace indigox::algorithm
