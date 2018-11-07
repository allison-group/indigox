#include <EASTL/vector_map.h>

#include "../../utils/fwd_declares.hpp"
#include "../../classes/parameterised.hpp"
#include "../../classes/athenaeum.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP
#define INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP

namespace indigox::algorithm {
  
  template <class V, class E, class S, class D, class VP, class EP>
  struct MappingCallback {
    using GraphType = graph::BaseGraph<V,E,S,D,VP,EP>;
    using CorrespondenceMap = eastl::vector_map<V, V>;
    virtual bool operator()(const CorrespondenceMap&) { return true; }
    virtual bool operator()(const V&, const V&) { return true; }
    virtual bool operator()(const E&, const E&) { return true; }
    virtual ~MappingCallback() = default;
  };
  
  struct CMGCallback
  : public MappingCallback<graph::CMGVertex, graph::CMGEdge, graph::sCondensedMolecularGraph,
  graph::Undirected, graph::GraphLabel, graph::GraphLabel>
  { };
  
  struct MGCallback
  : public MappingCallback<graph::MGVertex, graph::MGEdge, graph::sMolecularGraph,
  graph::Undirected, graph::GraphLabel, graph::GraphLabel>
  { };
  
  struct CMGPrintCallback : public CMGCallback {
    int count;
    CMGPrintCallback() : count(0) { }
    bool operator()(const CorrespondenceMap& cmap) override;
    bool operator()(const graph::CMGVertex& vs, const graph::CMGVertex& vl) override;
  };
  
  struct MGPrintCallback : public MGCallback {
    int count;
    MGPrintCallback() : count(0) { }
    bool operator()(const CorrespondenceMap& cmap) override {
      std::cout << "Mapping instance " << ++count << ":\n";
      for (auto& ab : cmap) std::cout << ab.first << " -- " << ab.second << "\n";
      std::cout << "\n";
      return true;
    }
  };

  void SubgraphIsomorphisms(graph::CondensedMolecularGraph& G1,
                            graph::CondensedMolecularGraph& G2,
                            CMGCallback& callback);
  
  void SubgraphIsomorphisms(graph::MolecularGraph& G1, graph::MolecularGraph& G2,
                            MGCallback& callback);
}


#endif /* INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP */
