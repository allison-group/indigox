#include <vector>

#include "../../graph/condensed.hpp"
#include "../../utils/common.hpp"
#include "../../utils/numerics.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP
#define INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP

namespace indigox::algorithm {
  
  /*! \brief Generate connected subgraphs in G.
   *  \param G the graph to generate connected subgraphs from.
   *  \param[out] subGs a vector to store the generated subgraphs in.
   *  \param minsize minimum size of subgraph to generate.
   *  \param maxsize maximum size of subgraph to generate.
   *  \return the number of generated subgraphs. */
  size_ ConnectedSubgraphs(const graph::CondensedMolecularGraph& G,
                           std::vector<graph::CondensedMolecularGraph>& subGs,
                           size_ minsize = 0, size_ maxsize = 1e18);
  
}

#endif /* INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP */
