#ifndef INDIGOX_ALGORITHM_ACCESS_HPP
#define INDIGOX_ALGORITHM_ACCESS_HPP

#include "../utils/fwd_declares.hpp"

#include "../graph/base_graph.hpp"

namespace indigox::algorithm {
  struct access {
    template <class V, class E, class S, class D, class VP, class EP>
    static typename graph::BaseGraph<V,E,S,D,VP,EP>::graph_type&
    GetGraph(graph::BaseGraph<V,E,S,D,VP,EP>& G) { return G._g; }
    
    template <class V, class E, class S, class D, class VP, class EP>
    static typename graph::BaseGraph<V,E,S,D,VP,EP>::VertMap&
    GetVertexMap(graph::BaseGraph<V,E,S,D,VP,EP>& G) { return G._vm; }
    
    template <class V, class E, class S, class D, class VP, class EP>
    static typename graph::BaseGraph<V,E,S,D,VP,EP>::EdgeMap&
    GetEdgeMap(graph::BaseGraph<V,E,S,D,VP,EP>& G) { return G._em; }
    
  };
}

#endif
