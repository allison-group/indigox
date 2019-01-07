#include "../algorithm/graph/connectivity.hpp"
#include "../algorithm/graph/cycles.hpp"
#include "../utils/serialise.hpp"
#include "../utils/triple.hpp"

#include <EASTL/vector_map.h>
#include <EASTL/vector_set.h>
#include <vector>

namespace indigox::graph {

#define GRAPHTEMPLATE                                                          \
  template <class V, class E, class S, class D, class VP, class EP>
#define BG BaseGraph<V, E, S, D, VP, EP>
#define BASEGRAPH(...) GRAPHTEMPLATE __VA_ARGS__ BG
#define tBG typename BG

  BASEGRAPH(struct)::BaseImpl {
    graph_type boost_graph;
    VertMap vertex_descriptors;
    EdgeMap edge_descriptors;
    VertContain vertices;
    EdgeContain edges;
    NbrsContain predecessors;
    NbrsContain successors;
    State state;

    ComponentContain cached_connected_components;
    State state_cached_components;

    VertContain cached_cyclic_vertices;
    EdgeContain cached_cyclic_edges;
    CycleEdgeContain cached_cycles;
    State state_cached_cycles;

    BaseImpl()
        : boost_graph(), state(0), state_cached_components(0),
          state_cached_cycles(0) {
    }

    V GetSourceVertex(const E &e) const {
      EdgeType eboost = edge_descriptors.left.at(e);
      VertType source = boost::source(eboost, boost_graph);
      return vertex_descriptors.right.at(source);
    }

    V GetTargetVertex(const E &e) const {
      EdgeType eboost = edge_descriptors.left.at(e);
      VertType target = boost::target(eboost, boost_graph);
      return vertex_descriptors.right.at(target);
    }

    void AddVertex(const V &v) {
      VertType vboost = boost::add_vertex(VP(), boost_graph);
      vertex_descriptors.insert(v, vboost);
      vertices.emplace_back(v);
      predecessors.emplace(v, VertContain());
      if (D::is_directed) {
        successors.emplace(v, VertContain());
      }
    }

    void AddEdge(const V &u, const V &v, const E &e) {
      // assume has both the vertices already. Make the BaseGraph method do the
      // real checking of this
      VertType uboost = vertex_descriptors.left.at(u);
      VertType vboost = vertex_descriptors.left.at(v);
      EdgeType eboost =
          boost::add_edge(uboost, vboost, EP(), boost_graph).first;
      edge_descriptors.insert(e, eboost);
      edges.emplace_back(e);
    }

    template <typename Archive>
    void save(Archive & archive, const uint32_t) const {
      std::vector<std::pair<V, V>> vertex_edges;
      vertex_edges.reserve(edges.size());
      for (const E &e : edges) {
        vertex_edges.emplace_back(GetSourceVertex(e), GetTargetVertex(e));
      }

      archive(INDIGOX_SERIAL_NVP("vertices", vertices),
              INDIGOX_SERIAL_NVP("vertex_edges", vertex_edges),
              INDIGOX_SERIAL_NVP("edges", edges),
              INDIGOX_SERIAL_NVP("state", state));
    }

    template <typename Archive> void load(Archive & archive, const uint32_t) {
      VertContain input_vertices;
      std::vector<std::pair<V, V>> input_vertex_edges;
      EdgeContain input_edges;
      archive(INDIGOX_SERIAL_NVP("vertices", input_vertices),
              INDIGOX_SERIAL_NVP("vertex_edges", input_vertex_edges),
              INDIGOX_SERIAL_NVP("edges", input_edges));

      // Build the graph
      for (V &v : input_vertices) {
        AddVertex(v);
      }
      for (size_t pos = 0; pos < input_edges.size(); ++pos) {
        AddEdge(input_vertex_edges[pos].first, input_vertex_edges[pos].second,
                input_edges[pos]);
      }

      // Reload state
      archive(INDIGOX_SERIAL_NVP("state", state));
    }
  };

  BASEGRAPH(template <typename Archive> void)::serialise(Archive &archive,
                                                         const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("base_graph_data", m_basedata));
  }

  //  BASEGRAPH(void)::Clear() {
  //    _g.clear();
  //    _vm.clear();
  //    _em.clear();
  //    _v.clear();
  //    _e.clear();
  //    _pre.clear();
  //    _suc.clear();
  //    _comp_cache.clear();
  //    _vcyclic_cache.clear();
  //    _ecyclic_cache.clear();
  //    _cycles_cache.clear();
  //  }

  BASEGRAPH(void)::AddVertex(const V &v) {
    ++m_basedata->state;
    m_basedata->AddVertex(v);
  }

  BASEGRAPH(void)::RemoveVertex(const V &v) {
    ++m_basedata->state;
    VertType vboost = GetDescriptor(v);
    // Remove adjacent edges
    PredIter vi, vi_end;
    std::tie(vi, vi_end) =
        boost::inv_adjacent_vertices(vboost, m_basedata->boost_graph);
    for (; vi != vi_end; ++vi) {
      V u = GetV(*vi);
      E e = GetE(boost::edge(vboost, *vi, m_basedata->boost_graph).first);
      m_basedata->edge_descriptors.erase(e);
      m_basedata->edges.erase(
          std::find(m_basedata->edges.begin(), m_basedata->edges.end(), e));
    }
    // Remove incident edges of directed graphs
    if (D::is_directed) {
      NbrsIter vp, vp_end;
      std::tie(vp, vp_end) =
          boost::adjacent_vertices(vboost, m_basedata->boost_graph);
      for (; vp != vp_end; ++vp) {
        V u = GetV(*vp);
        E e = GetE(boost::edge(*vp, vboost, m_basedata->boost_graph).first);
        m_basedata->edge_descriptors.erase(e);
        m_basedata->edges.erase(
            std::find(m_basedata->edges.begin(), m_basedata->edges.end(), e));
      }
    }
    // Remove the vertex
    m_basedata->vertex_descriptors.erase(v);
    m_basedata->vertices.erase(
        std::find(m_basedata->vertices.begin(), m_basedata->vertices.end(), v));
    m_basedata->predecessors.erase(v);
    if (D::is_directed)
      m_basedata->successors.erase(v);
    boost::clear_vertex(vboost, m_basedata->boost_graph);
    boost::remove_vertex(vboost, m_basedata->boost_graph);
  }

  BASEGRAPH(void)::AddEdge(const V &u, const V &v, const E &e) {
    ++m_basedata->state;
    if (!HasVertex(u))
      AddVertex(u);
    if (!HasVertex(v))
      AddVertex(v);
    m_basedata->AddEdge(u, v, e);
  }

  BASEGRAPH(void)::RemoveEdge(const E &e) {
    ++m_basedata->state;
    EdgeType eboost = GetDescriptor(e);
    m_basedata->edge_descriptors.erase(eboost);
    m_basedata->edges.erase(
        std::find(m_basedata->edges.begin(), m_basedata->edges.end(), e));
    boost::remove_edge(eboost, m_basedata->boost_graph);
  }

  BASEGRAPH(void)::RemoveEdge(const V &u, const V &v) {
    E e = GetEdge(u, v);
    RemoveEdge(e);
  }

  BASEGRAPH(bool)::HasVertex(const V &v) const {
    return m_basedata->vertex_descriptors.left.find(v) !=
           m_basedata->vertex_descriptors.left.end();
  }

  BASEGRAPH(bool)::HasEdge(const E &e) const {
    return m_basedata->edge_descriptors.left.find(e) !=
           m_basedata->edge_descriptors.left.end();
  }

  BASEGRAPH(bool)::HasEdge(const V &u, const V &v) const {
    if (!HasVertex(u) || !HasVertex(v))
      return false;
    VertType u_ = GetDescriptor(u);
    VertType v_ = GetDescriptor(v);
    return boost::edge(u_, v_, m_basedata->boost_graph).second;
  }

  BASEGRAPH(int64_t)::NumVertices() const {
    return boost::num_vertices(m_basedata->boost_graph);
  }

  BASEGRAPH(int64_t)::NumEdges() const {
    return boost::num_edges(m_basedata->boost_graph);
  }

  BASEGRAPH(int64_t)::Degree(const V &v) const {
    return HasVertex(v) ? OutDegree(GetDescriptor(v)) : -1;
  }

  BASEGRAPH(int64_t)::InDegree(const V &v) const {
    if (!HasVertex(v))
      return -1;
    if (D::is_directed)
      return InDegree(GetDescriptor(v));
    return OutDegree(GetDescriptor(v));
  }

  BASEGRAPH(const tBG::VertContain &)::GetNeighbours(const V &v) {
    static_assert(!D::is_directed, "Requires an undirected graph.");
    VertType vboost = GetDescriptor(v);
    auto adjis = boost::adjacent_vertices(vboost, m_basedata->boost_graph);
    m_basedata->predecessors.at(v).clear();
    for (; adjis.first != adjis.second; ++adjis.first)
      m_basedata->predecessors.at(v).emplace_back(GetV(*adjis.first));
    return m_basedata->predecessors.at(v);
  }

  BASEGRAPH(const tBG::VertContain &)::GetPredecessors(const V &v) {
    static_assert(D::is_directed, "Requires a directed graph.");
    VertType vboost = GetDescriptor(v);
    auto adjis = boost::inv_adjacent_vertices(vboost, m_basedata->boost_graph);
    m_basedata->predecessors.at(v).clear();
    for (; adjis.first != adjis.second; ++adjis.first)
      m_basedata->predecessors.at(v).emplace_back(GetV(*adjis.first));
    return m_basedata->predecessors.at(v);
  }

  BASEGRAPH(const tBG::VertContain &)::GetSuccessors(const V &v) {
    static_assert(D::is_directed, "Requires a directed graph.");
    VertType vboost = GetDescriptor(v);
    auto adjis = boost::adjacent_vertices(vboost, m_basedata->boost_graph);
    m_basedata->successors.at(v).clear();
    for (; adjis.first != adjis.second; ++adjis.first)
      m_basedata->successors.at(v).emplace_back(GetV(*adjis.first));
    return m_basedata->successors.at(v);
  }

  BASEGRAPH(std::pair<V, V>)::GetVertices(const E &e) const {
    EdgeType eboost = GetDescriptor(e);
    VertType u = boost::source(eboost, m_basedata->boost_graph);
    VertType v = boost::target(eboost, m_basedata->boost_graph);
    return {GetV(u), GetV(v)};
  }

  BASEGRAPH(const tBG::VertContain &)::GetVertices() const {
    return m_basedata->vertices;
  }

  BASEGRAPH(const tBG::EdgeContain &)::GetEdges() const {
    return m_basedata->edges;
  }

  BASEGRAPH(E)::GetEdge(const V &u, const V &v) const {
    VertType uboost = GetDescriptor(u);
    VertType vboost = GetDescriptor(v);
    EdgeType eboost =
        boost::edge(uboost, vboost, m_basedata->boost_graph).first;
    return GetE(eboost);
  }

  BASEGRAPH(V)::GetSourceVertex(const E &e) const {
    return m_basedata->GetSourceVertex(e);
  }

  BASEGRAPH(V)::GetTargetVertex(const E &e) const {
    return m_basedata->GetTargetVertex(e);
  }

  // =====================================================================
  // == Private helper functions =========================================
  // =====================================================================
  BASEGRAPH(tBG::VertType)::GetDescriptor(const V &v) const {
    return m_basedata->vertex_descriptors.left.at(v);
  }

  BASEGRAPH(V)::GetV(BG::VertType v) const {
    return m_basedata->vertex_descriptors.right.at(v);
  }

  BASEGRAPH(tBG::EdgeType)::GetDescriptor(const E &e) const {
    return m_basedata->edge_descriptors.left.at(e);
  }

  BASEGRAPH(E)::GetE(BG::EdgeType e) const {
    return m_basedata->edge_descriptors.right.at(e);
  }

  BASEGRAPH(int64_t)::OutDegree(BG::VertType v) const {
    return boost::out_degree(v, m_basedata->boost_graph);
  }

  BASEGRAPH(int64_t)::InDegree(BG::VertType v) const {
    return boost::in_degree(v, m_basedata->boost_graph);
  }

  // ======================================================================
  // == Algorithmic stuffs ================================================
  // ======================================================================

  // Connectivity
  BASEGRAPH(bool)::IsConnected() {
    GetConnectedComponents();
    return m_basedata->cached_connected_components.size() == 1;
  }

  BASEGRAPH(int64_t)::NumConnectedComponents() {
    GetConnectedComponents();
    return m_basedata->cached_connected_components.size();
  }

  BASEGRAPH(const tBG::ComponentContain &)::GetConnectedComponents() {
    if (m_basedata->state_cached_components == m_basedata->state)
      return m_basedata->cached_connected_components;
    algorithm::ConnectedComponents(*this,
                                   m_basedata->cached_connected_components);
    m_basedata->state_cached_components = m_basedata->state;
    return m_basedata->cached_connected_components;
  }

  // Cycles
  BASEGRAPH(bool)::IsCyclic(const V &v) {
    GetCycles();
    auto pos = std::find(m_basedata->cached_cyclic_vertices.begin(),
                         m_basedata->cached_cyclic_vertices.end(), v);
    return pos != m_basedata->cached_cyclic_vertices.end();
  }

  BASEGRAPH(bool)::IsCyclic(const V &v, uint32_t sz) {
    if (!IsCyclic(v))
      return false;
    for (auto &cyc : m_basedata->cached_cycles) {
      if (cyc.size() > sz)
        return false;
      for (E &edge : cyc) {
        V a = GetSourceVertex(edge);
        V b = GetTargetVertex(edge);
        if (a == v || b == v)
          return true;
      }
    }
    return false;
  }

  BASEGRAPH(bool)::IsCyclic(const E &e) {
    GetCycles();
    auto pos = std::find(m_basedata->cached_cyclic_edges.begin(),
                         m_basedata->cached_cyclic_edges.end(), e);
    return pos != m_basedata->cached_cyclic_edges.end();
  }

  BASEGRAPH(bool)::IsCyclic(const E &e, uint32_t sz) {
    if (!IsCyclic(e))
      return false;
    for (auto &cyc : m_basedata->cached_cycles) {
      if (cyc.size() > sz)
        return false;
      for (E &edge : cyc) {
        if (edge == e)
          return true;
      }
    }
    return false;
  }

  BASEGRAPH(const tBG::CycleEdgeContain &)::GetCycles() {
    using VertSet = eastl::vector_set<V>;
    using EdgeSet = eastl::vector_set<E>;

    if (m_basedata->state == m_basedata->state_cached_cycles)
      return m_basedata->cached_cycles;
    algorithm::AllCycles(*this, m_basedata->cached_cycles);
    VertSet cyclic_v;
    cyclic_v.reserve(NumVertices());
    EdgeSet cyclic_e;
    cyclic_e.reserve(NumEdges());

    for (EdgeContain &cycle : m_basedata->cached_cycles) {
      for (E &edge : cycle) {
        cyclic_e.emplace(edge);
        cyclic_v.emplace(GetSourceVertex(edge));
        cyclic_v.emplace(GetTargetVertex(edge));
      }
    }
    m_basedata->cached_cyclic_vertices.clear();
    m_basedata->cached_cyclic_vertices.assign(cyclic_v.begin(), cyclic_v.end());
    m_basedata->cached_cyclic_edges.clear();
    m_basedata->cached_cyclic_edges.assign(cyclic_e.begin(), cyclic_e.end());
    std::sort(
        m_basedata->cached_cycles.begin(), m_basedata->cached_cycles.end(),
        [](EdgeContain &a, EdgeContain &b) { return a.size() < b.size(); });
    m_basedata->state_cached_cycles = m_basedata->state;
    return m_basedata->cached_cycles;
  }

  BASEGRAPH(int64_t)::NumCycles() {
    GetCycles();
    return m_basedata->cached_cycles.size();
  }

  // Isomorphism

#undef GRAPHTEMPLATE
#undef BASEGRAPH
#undef BG
#undef tBG
} // namespace indigox::graph
