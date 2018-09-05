#include <map>
#include <vector>

#include "../../utils/fwd_declares.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP
#define INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP

namespace indigox::algorithm {
  
  void SubgraphIsomorphisms(graph::CondensedMolecularGraph G1,
                            graph::CondensedMolecularGraph G2,
                            std::map<graph::CMGVertex,
                            std::vector<graph::CMGVertex>>& maps);
  
  void SubgraphIsomorphisms(graph::MolecularGraph G1,
                            graph::MolecularGraph G2,
                            std::map<graph::MGVertex,
                            std::vector<graph::MGVertex>>& maps);
  
}

#endif /* INDIGOX_ALGORITHM_GRAPH_ISOMORPHISM_HPP */
