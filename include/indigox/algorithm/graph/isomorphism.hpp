#include "../../utils/fwd_declares.hpp"

#include <EASTL/vector_map.h>
#include <EASTL/bitset.h>
#include <rilib/RI.h>

#ifndef INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP
#define INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP

namespace indigox::graph {
  //! \brief type used to store the isomorphism testing mask for IXCMGVertex.
  using VertexIsoMask = eastl::bitset<37, uint64_t>;
  //! \brief type used to store the isomorphism testing mask for IXCMGEdge.
  using EdgeIsoMask = eastl::bitset<14, uint16_t>;
}

namespace indigox::algorithm {

  template <class V, class E, class S, class D, class VP, class EP>
  struct MappingCallback {
    using GraphType = graph::BaseGraph<V, E, S, D, VP, EP>;
    using CorrespondenceMap = eastl::vector_map<V, V>;
    virtual bool operator()(const CorrespondenceMap &) { return true; }
    virtual bool operator()(const V &, const V &) { return true; }
    virtual bool operator()(const E &, const E &) { return true; }
    virtual ~MappingCallback() = default;
  };

  struct CMGCallback
      : public MappingCallback<
            graph::CMGVertex, graph::CMGEdge, graph::CondensedMolecularGraph,
            graph::Undirected, graph::GraphLabel, graph::GraphLabel> {};

  struct MGCallback
      : public MappingCallback<graph::MGVertex, graph::MGEdge,
                               graph::MolecularGraph, graph::Undirected,
                               graph::GraphLabel, graph::GraphLabel> {};

  struct CMGPrintCallback : public CMGCallback {
    int count;
    CMGPrintCallback() : count(0) {}
    bool operator()(const CorrespondenceMap &cmap) override;
    bool operator()(const graph::CMGVertex &vs,
                    const graph::CMGVertex &vl) override;
  };

  struct MGPrintCallback : public MGCallback {
    int count;
    MGPrintCallback() : count(0) {}
    bool operator()(const CorrespondenceMap &cmap) override;
  };

  void SubgraphIsomorphisms(graph::CondensedMolecularGraph &G1,
                            graph::CondensedMolecularGraph &G2,
                            CMGCallback &callback);

  void SubgraphIsomorphisms(graph::MolecularGraph &G1,
                            graph::MolecularGraph &G2, MGCallback &callback);
  
  int LargestCommonSubgraph(graph::CondensedMolecularGraph& source_g,
                            graph::CondensedMolecularGraph& target_g,
                            size_t smallest_size,
                            eastl::vector_map<graph::CMGVertex, eastl::vector<std::pair<graph::CMGVertex, graph::CMGVertex>>>& per_vertex_largest_subgraphs);
  
  int LargestCommonSubgraph(graph::MolecularGraph& source_g,
                            graph::MolecularGraph& target_g,
                            size_t smallest_size,
                            eastl::vector_map<graph::MGVertex, eastl::vector<std::pair<graph::MGVertex, graph::MGVertex>>>& per_vertex_largest_subgraphs);
  
  struct Uint64AttrComparator : public rilib::AttributeComparator {
    Uint64AttrComparator(){};
    virtual bool compare(void *attr1, void *attr2) {
      return (*((uint64_t *)attr1) - *((uint64_t *)attr2)) == 0;
    }
    virtual int compareint(void *attr1, void *attr2) {
      return *((uint64_t *)attr1) - *((uint64_t *)attr2);
    }
  };
  
  struct Uint32AttrComparator : public rilib::AttributeComparator {
    Uint32AttrComparator(){};
    virtual bool compare(void *attr1, void *attr2) {
      return (*((uint32_t *)attr1) - *((uint32_t *)attr2)) == 0;
    }
    virtual int compareint(void *attr1, void *attr2) {
      return *((uint32_t *)attr1) - *((uint32_t *)attr2);
    }
  };
  
  std::unique_ptr<rilib::Graph>
  CMGToRIGraph(graph::CondensedMolecularGraph &cmg, graph::EdgeIsoMask emask,
               graph::VertexIsoMask vmask);
  std::unique_ptr<rilib::Graph> MGToRIGraph(graph::MolecularGraph& mg);
  
} // namespace indigox::algorithm

#endif /* INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP */
