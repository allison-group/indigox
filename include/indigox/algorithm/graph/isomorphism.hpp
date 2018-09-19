#include <map>
#include <vector>

#include <EASTL/vector_map.h>

#include "../../utils/fwd_declares.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP
#define INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP

namespace indigox::algorithm {
  
  template <class GraphType>
  struct MappingCallback {
    using Vertex = typename GraphType::VertexType;
    using CorrespondenceMap = eastl::vector_map<Vertex, Vertex>;
    
    virtual bool operator()(const CorrespondenceMap&) { return true; }
    
  };
  
  void SubgraphIsomorphisms(graph::CondensedMolecularGraph G1,
                            graph::CondensedMolecularGraph G2,
                            MappingCallback<graph::IXCondensedMolecularGraph>& callback);
  
  void SubgraphIsomorphisms(graph::MolecularGraph G1,
                            graph::MolecularGraph G2,
                            MappingCallback<graph::IXMolecularGraph>& callback);
  
}

#endif /* INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP */
