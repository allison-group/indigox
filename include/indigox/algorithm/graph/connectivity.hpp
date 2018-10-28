#ifndef INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP
#define INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP

#include <cstdint>
#include <vector>

#include <boost/graph/connected_components.hpp>

#include "../access.hpp"
#include "../../utils/common.hpp"
#include "../../utils/fwd_declares.hpp"

#include "../../graph/base_graph.hpp"

namespace indigox::algorithm {
  
  template<class V, class E, class D, class VP, class EP>
  int64_t ConnectedComponents(graph::BaseGraph<V,E,D,VP,EP>& G,
                              typename graph::BaseGraph<V,E,D,VP,EP>::ComponentContain& bits) {
    static_assert(!D::is_directed, "Connectivity requires an undirected graph");
    
    using GraphType = graph::BaseGraph<V,E,D,VP,EP>;
    using BaseType = typename GraphType::graph_type;
    using VertIdxMap = eastl::vector_map<typename GraphType::VertType, int>;
    using namespace boost;
    
    BaseType& g = access::GetGraph(G);
    VertIdxMap data;
    associative_property_map<VertIdxMap> indexMap(data);
    typename GraphType::VertIter begin, end;
    std::tie(begin, end) = boost::vertices(g);
    for (int i = 0; begin != end; ++begin, ++i) put(indexMap, *begin, i);
    size_t num = connected_components(g, get(&VP::component, g),
                                      vertex_index_map(indexMap));
    bits.clear();
    bits.assign(num, typename GraphType::VertContain());
    for (auto& v : access::GetVertexMap(G).right)
      bits[g[v.first].component].emplace_back(v.second);
    return bits.size();
  }
  
  /*! \brief Generate connected subgraphs in G.
   *  \param G the graph to generate connected subgraphs from.
   *  \param[out] subGs a vector to store the generated subgraphs in.
   *  \param minsize minimum size of subgraph to generate.
   *  \param maxsize maximum size of subgraph to generate.
   *  \return the number of generated subgraphs. */
  size_t ConnectedSubgraphs(const graph::CondensedMolecularGraph& G,
                           std::vector<graph::CondensedMolecularGraph>& subGs,
                           size_t minsize = 0, size_t maxsize = 1e18);
  
  /*! \brief Generate connected subgraphs in G.
   *  \param G the graph to generate connected subgraphs from.
   *  \param[out] subGs a vector to store the generated subgraphs in.
   *  \param minsize minimum size of subgraph to generate.
   *  \param maxsize maximum size of subgraph to generate.
   *  \return the number of generated subgraphs. */
  size_t ConnectedSubgraphs(const graph::MolecularGraph& G,
                           std::vector<graph::MolecularGraph>& subGs,
                           size_t minsize = 0, size_t maxsize = 1e18);
  
}

#endif /* INDIGOX_ALGORITHM_GRAPH_CONNECTIVITY_HPP */
