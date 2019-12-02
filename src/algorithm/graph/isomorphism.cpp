#include <indigox/algorithm/access.hpp>
#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/classes/periodictable.hpp>

#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <rilib/RI.h>

#include <EASTL/vector_map.h>


namespace indigox::algorithm {
  using namespace indigox::graph;
  
  std::unique_ptr<rilib::Graph>
  CMGToRIGraph(graph::CondensedMolecularGraph &cmg, graph::EdgeIsoMask emask,
               graph::VertexIsoMask vmask) {
    std::unique_ptr<rilib::Graph> G = std::make_unique<rilib::Graph>();
    
    std::vector<graph::CMGVertex> cmg_v = cmg.GetVertices();
    
    // Build the vertices
    G->nof_nodes = cmg_v.size();
    G->nodes_attrs = (void **)malloc(G->nof_nodes * sizeof(void *));
    G->out_adj_sizes = (int *)calloc(G->nof_nodes, sizeof(int));
    G->in_adj_sizes = (int *)calloc(G->nof_nodes, sizeof(int));
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->nodes_attrs[i] = (uint64_t *)malloc(sizeof(uint64_t));
      *((uint64_t *)G->nodes_attrs[i]) =
      (cmg_v[i].GetIsomorphismMask() & vmask).to_uint64();
      G->out_adj_sizes[i] += cmg.Degree(cmg_v[i]);
      G->in_adj_sizes[i] += cmg.Degree(cmg_v[i]);
    }
    
    // Build the edges
    G->out_adj_list = (int **)malloc(G->nof_nodes * sizeof(int *));
    G->in_adj_list = (int **)malloc(G->nof_nodes * sizeof(int *));
    G->out_adj_attrs = (void ***)malloc(G->nof_nodes * sizeof(void **));
    
    int *ink = (int *)calloc(G->nof_nodes, sizeof(int));
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->in_adj_list[i] = (int *)calloc(G->in_adj_sizes[i], sizeof(int));
    }
    
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->out_adj_list[i] = (int *)calloc(G->in_adj_sizes[i], sizeof(int));
      G->out_adj_attrs[i] =
      (void **)malloc(G->out_adj_sizes[i] * sizeof(void *));
      std::vector<graph::CMGVertex> nbrs = cmg.GetNeighbours(cmg_v[i]);
      for (int j = 0; j < G->out_adj_sizes[i]; ++j) {
        int idx = std::distance(cmg_v.begin(),
                                std::find(cmg_v.begin(), cmg_v.end(), nbrs[j]));
        G->out_adj_list[i][j] = idx;
        G->out_adj_attrs[i][j] = (uint32_t *)malloc(sizeof(uint32_t));
        *((uint32_t *)G->out_adj_attrs[i][j]) =
        (cmg.GetEdge(cmg_v[i], nbrs[j]).GetIsomorphismMask() & emask)
        .to_uint32();
        G->in_adj_list[idx][ink[idx]] = i;
        ink[idx]++;
      }
    }
    
    free(ink);
    return G;
  }
  
  std::unique_ptr<rilib::Graph> MGToRIGraph(graph::MolecularGraph& mg) {
    std::unique_ptr<rilib::Graph> G = std::make_unique<rilib::Graph>();
    std::vector<graph::MGVertex> mg_v = mg.GetVertices();
    
    // Build the vertices
    G->nof_nodes = mg_v.size();
    G->nodes_attrs = (void **)malloc(G->nof_nodes * sizeof(void *));
    G->out_adj_sizes = (int *)calloc(G->nof_nodes, sizeof(int));
    G->in_adj_sizes = (int *)calloc(G->nof_nodes, sizeof(int));
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->nodes_attrs[i] = (uint32_t *)malloc(sizeof(uint32_t));
      *((uint32_t *)G->nodes_attrs[i]) = mg_v[i].GetAtom().GetElement().GetAtomicNumber();
      G->out_adj_sizes[i] += mg.Degree(mg_v[i]);
      G->in_adj_sizes[i] += mg.Degree(mg_v[i]);
    }
    
    // Build the edges
    G->out_adj_list = (int **)malloc(G->nof_nodes * sizeof(int *));
    G->in_adj_list = (int **)malloc(G->nof_nodes * sizeof(int *));
    G->out_adj_attrs = (void ***)malloc(G->nof_nodes * sizeof(void **));
    
    int *ink = (int *)calloc(G->nof_nodes, sizeof(int));
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->in_adj_list[i] = (int *)calloc(G->in_adj_sizes[i], sizeof(int));
    }
    
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->out_adj_list[i] = (int *)calloc(G->in_adj_sizes[i], sizeof(int));
      G->out_adj_attrs[i] = (void **)malloc(G->out_adj_sizes[i] * sizeof(void *));
      std::vector<graph::MGVertex> nbrs = mg.GetNeighbours(mg_v[i]);
      for (int j = 0; j < G->out_adj_sizes[i]; ++j) {
        int idx = std::distance(mg_v.begin(),
                                std::find(mg_v.begin(), mg_v.end(), nbrs[j]));
        G->out_adj_list[i][j] = idx;
        G->out_adj_attrs[i][j] = (uint32_t *)malloc(sizeof(uint32_t));
        *((uint32_t *)G->out_adj_attrs[i][j]) = (uint32_t)mg.GetEdge(mg_v[i], nbrs[j]).GetBond().GetOrder();
        G->in_adj_list[idx][ink[idx]] = i;
        ink[idx]++;
      }
    }
    
    free(ink);
    
    return G;
  }
  
  // =========================================================================
  // ==== Boost VF2 based subgraph isomorphism implementation ================
  // =========================================================================

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
  void VF2SubgraphIsoRunner(graph::BaseGraph<V, E, S, D, VP, EP> &G1,
                            graph::BaseGraph<V, E, S, D, VP, EP> &G2,
                            MappingCallback<V, E, S, D, VP, EP> &callback) {

    using GraphType = graph::BaseGraph<V, E, S, D, VP, EP>;
    using BoostGraph = typename GraphType::graph_type;
    using Vertex = typename GraphType::VertType;
    using VertIter = typename GraphType::VertIter;
    using VertexMap = eastl::vector_map<Vertex, size_t>;

    if (G2.NumVertices() < G1.NumVertices()) {
      VF2SubgraphIsoRunner(G2, G1, callback);
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
    VF2SubgraphIsoRunner(G1, G2, CB);
  }

  void SubgraphIsomorphisms(MolecularGraph &G1, MolecularGraph &G2,
                            MGCallback &CB) {
    VF2SubgraphIsoRunner(G1, G2, CB);
  }

  using CSubgraphMap = eastl::vector<std::pair<CMGVertex, CMGVertex>>;
  using AllCSubgraphMaps = eastl::vector_map<CMGVertex, CSubgraphMap>;
  struct RILargestCommonCSubgraphMatcher : rilib::MatchListener {
    CondensedMolecularGraph &small, &big;
    AllCSubgraphMaps& matches;
    RILargestCommonCSubgraphMatcher(CondensedMolecularGraph& s, CondensedMolecularGraph& b,
                                   AllCSubgraphMaps& m)
    : small(s), big(b), matches(m) { }
    virtual void match(int n, int * small_ids, int * large_ids) {
      CSubgraphMap m; m.reserve(n);
      for (int i = 0; i < n; ++i) {
        m.emplace_back(small.GetVertices()[small_ids[i]],
                                    big.GetVertices()[large_ids[i]]);
      }
      for (auto& sub_tar : m) {
        if (matches.at(sub_tar.first).size() < m.size())
          matches[sub_tar.first] = m;
      }
    }
  };
  
  // LCSubGraph
  int LargestCommonSubgraph(CondensedMolecularGraph& source_g,
                            CondensedMolecularGraph& target_g,
                            size_t smallest_size,
                            AllCSubgraphMaps& largest_subgraphs) {
    
    // Setup the return map
    largest_subgraphs.clear();
    for(CMGVertex v : source_g.GetVertices()) largest_subgraphs[v] = CSubgraphMap();
    
    // Generate all subgraphs of source_g
    ConnectedSubgraphs subgraph_generator(source_g, smallest_size, std::min(target_g.NumVertices(), source_g.NumVertices()));
    CondensedMolecularGraph sub;
    graph::VertexIsoMask vertmask; vertmask.set();
    graph::EdgeIsoMask edgemask; edgemask.set();
    std::unique_ptr<rilib::Graph> target_ri = CMGToRIGraph(target_g, edgemask, vertmask);
    // Iterate over subgraphs checking for isomporhisms
    while (subgraph_generator(sub)) {
      // convert subgraph to ri graph
      std::unique_ptr<rilib::Graph> sub_ri = CMGToRIGraph(sub, edgemask, vertmask);
      Uint64AttrComparator* vert_compare = new Uint64AttrComparator();
      Uint32AttrComparator* edge_compare = new Uint32AttrComparator();
      RILargestCommonCSubgraphMatcher* listener = new RILargestCommonCSubgraphMatcher(sub, target_g, largest_subgraphs);
      rilib::MaMaConstrFirst* mama = new rilib::MaMaConstrFirst(*sub_ri);
      mama->build(*sub_ri);
      long tmp_1, tmp_2, tmp_3;
      // run the matching
      rilib::match(*target_ri, *sub_ri, *mama, *listener,
                   rilib::MATCH_TYPE::MT_INDSUB, *vert_compare,
                   *edge_compare, &tmp_1, &tmp_2, &tmp_3);
      delete vert_compare;
      delete edge_compare;
      delete listener;
      delete mama;
    }
    
    return 0;
  }
  
  using SubgraphMap = eastl::vector<std::pair<MGVertex, MGVertex>>;
  using AllSubgraphMaps = eastl::vector_map<MGVertex, SubgraphMap>;
  struct RILargestCommonSubgraphMatcher : rilib::MatchListener {
    MolecularGraph &small, &big;
    AllSubgraphMaps& matches;
    RILargestCommonSubgraphMatcher(MolecularGraph& s, MolecularGraph& b,
                                    AllSubgraphMaps& m)
    : small(s), big(b), matches(m) { }
    virtual void match(int n, int * small_ids, int * large_ids) {
      SubgraphMap m; m.reserve(n);
      for (int i = 0; i < n; ++i) {
        m.emplace_back(small.GetVertices()[small_ids[i]],
                       big.GetVertices()[large_ids[i]]);
      }
      for (auto& sub_tar : m) {
        if (matches.at(sub_tar.first).size() < m.size())
          matches[sub_tar.first] = m;
      }
    }
  };
  
  int LargestCommonSubgraph(MolecularGraph& source_g, MolecularGraph& target_g,
                            size_t smallest_size, AllSubgraphMaps& largest_subgraphs) {
    // Setup the return map
    largest_subgraphs.clear();
    for(MGVertex v : source_g.GetVertices()) largest_subgraphs[v] = SubgraphMap();
    
    // Generate all subgraphs of source_g
    ConnectedSubgraphs subgraph_generator(source_g, smallest_size, std::min(target_g.NumVertices(), source_g.NumVertices()));
    MolecularGraph sub;
    std::unique_ptr<rilib::Graph> target_ri = MGToRIGraph(target_g);
    // Iterate over subgraphs checking for isomporhisms
    while (subgraph_generator(sub)) {
      // convert subgraph to ri graph
      std::unique_ptr<rilib::Graph> sub_ri = MGToRIGraph(sub);
      Uint32AttrComparator* vert_compare = new Uint32AttrComparator();
      Uint32AttrComparator* edge_compare = new Uint32AttrComparator();
      RILargestCommonSubgraphMatcher* listener = new RILargestCommonSubgraphMatcher(sub, target_g, largest_subgraphs);
      rilib::MaMaConstrFirst* mama = new rilib::MaMaConstrFirst(*sub_ri);
      mama->build(*sub_ri);
      long tmp_1, tmp_2, tmp_3;
      // run the matching
      rilib::match(*target_ri, *sub_ri, *mama, *listener,
                   rilib::MATCH_TYPE::MT_INDSUB, *vert_compare,
                   *edge_compare, &tmp_1, &tmp_2, &tmp_3);
      delete vert_compare;
      delete edge_compare;
      delete listener;
      delete mama;
    }
    
    return 0;
  }
  
} // namespace indigox::algorithm
