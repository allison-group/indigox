#ifndef INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP
#define INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP

#include <cstdint>
#include <limits>
#include <vector>

#include "../../utils/fwd_declares.hpp"

namespace indigox::algorithm {
  
  template<class V, class E, class S, class D, class VP, class EP, class Container>
  int64_t ConnectedComponents(graph::BaseGraph<V,E,S,D,VP,EP>& G, Container& bits);
  
  /*! \brief Generate connected subgraphs in G.
   *  \param G the graph to generate connected subgraphs from.
   *  \param[out] subGs a vector to store the generated subgraphs in.
   *  \param minsize minimum size of subgraph to generate.
   *  \param maxsize maximum size of subgraph to generate.
   *  \return the number of generated subgraphs. */
  template<class V, class E, class S, class D, class VP, class EP>
  int64_t ConnectedSubgraphs(graph::BaseGraph<V,E,S,D,VP,EP>& G,
                             std::vector<S>& sub_graphs, size_t min = 0,
                             size_t max = std::numeric_limits<size_t>::max());
}

#endif /* INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP */
