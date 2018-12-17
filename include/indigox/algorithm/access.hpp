#ifndef INDIGOX_ALGORITHM_ACCESS_HPP
#define INDIGOX_ALGORITHM_ACCESS_HPP

#include "../graph/base_graph.hpp"
#include "../utils/fwd_declares.hpp"

namespace indigox::algorithm {
  struct access {
    template <class V, class E, class S, class D, class VP, class EP>
    static typename graph::BaseGraph<V, E, S, D, VP, EP>::graph_type &
    GetGraph(graph::BaseGraph<V, E, S, D, VP, EP> &G) {
      return G.m_basedata->boost_graph;
    }

    template <class V, class E, class S, class D, class VP, class EP>
    static typename graph::BaseGraph<V, E, S, D, VP, EP>::VertMap &
    GetVertexMap(graph::BaseGraph<V, E, S, D, VP, EP> &G) {
      return G.m_basedata->vertex_descriptors;
    }

    template <class V, class E, class S, class D, class VP, class EP>
    static typename graph::BaseGraph<V, E, S, D, VP, EP>::EdgeMap &
    GetEdgeMap(graph::BaseGraph<V, E, S, D, VP, EP> &G) {
      return G.m_basedata->edge_descriptors;
    }
  };
} // namespace indigox::algorithm

#endif
