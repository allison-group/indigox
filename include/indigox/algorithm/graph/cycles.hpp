#include "../../utils/fwd_declares.hpp"
#include "../../utils/numerics.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP
#define INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP

namespace indigox::algorithm {
  
  template<class V, class E, class S, class D, class VP, class EP, class Container>
  int64_t CycleBasis(graph::BaseGraph<V,E,S,D,VP,EP>& G, Container& basis);
  
  template<class V, class E, class S, class D, class VP, class EP, class Container>
  int64_t AllCycles(graph::BaseGraph<V,E,S,D,VP,EP>& G, Container& ecycles);
  
}

#endif /* INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP */
