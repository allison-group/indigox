#include <vector>

#include <EASTL/vector_map.h>
#include <EASTL/vector_set.h>

#include "../algorithm/graph/connectivity.hpp"
#include "../algorithm/graph/cycles.hpp"
#include "../utils/serialise.hpp"
#include "../utils/triple.hpp"

namespace indigox::graph {
  
#define GRAPHTEMPLATE template<class V, class E, class D, class VP, class EP>
#define BG BaseGraph<V,E,D,VP,EP>
#define BASEGRAPH(...) GRAPHTEMPLATE __VA_ARGS__ BG
#define tBG typename BG
  
  BASEGRAPH(template <typename Archive> void)::save(Archive& archive, const uint32_t) const {
    std::vector<stdx::triple<V, V, E>> edges;
    edges.reserve(NumEdges());
    for (const E& e : _e)
      edges.emplace_back(GetSourceVertex(e), GetTargetVertex(e), e);
    
    archive(INDIGOX_SERIAL_NVP("vertices", _v),
            INDIGOX_SERIAL_NVP("edges", edges),
            INDIGOX_SERIAL_NVP("state", cereal::base_class<utils::ModifiableObject>(this)));
  }
  
  BASEGRAPH(template <typename Archive> void)::load(Archive& archive, const uint32_t) {
    VertContain vertices;
    std::vector<stdx::triple<V, V, E>> edges;
    archive(INDIGOX_SERIAL_NVP("vertices", vertices),
            INDIGOX_SERIAL_NVP("edges", edges));
    
    // Build the graph
    for (V& v : vertices) AddVertex(v);
    for (auto& e : edges) AddEdge(e.first, e.second, e.third);
    
    // Reload state
    archive(INDIGOX_SERIAL_NVP("state", cereal::base_class<utils::ModifiableObject>(this)));
  }
  
  BASEGRAPH(void)::Clear() {
    _g.clear();
    _vm.clear();
    _em.clear();
    _v.clear();
    _e.clear();
    _pre.clear();
    _suc.clear();
    _comp_cache.clear();
    _vcyclic_cache.clear();
    _ecyclic_cache.clear();
    _cycles_cache.clear();
  }
  
  BASEGRAPH(void)::AddVertex(const V& v)  {
    ModificationMade();
    VertType vboost = boost::add_vertex(VP(), _g);
    _vm.insert(v, vboost);
    _v.emplace_back(v);
    _pre.emplace(v, VertContain());
    if (D::is_directed) _suc.emplace(v, VertContain());
  }
  
  BASEGRAPH(void)::RemoveVertex(const V& v) {
    ModificationMade();
    VertType vboost = GetDescriptor(v);
    // Remove adjacent edges
    PredIter vi, vi_end;
    std::tie(vi, vi_end) = boost::inv_adjacent_vertices(vboost, _g);
    for (; vi != vi_end; ++vi) {
      V u = GetV(*vi);
      E e = GetE(boost::edge(vboost, *vi, _g).first);
      _em.erase(e);
      _e.erase(std::find(_e.begin(), _e.end(), e));
    }
    // Remove incident edges of directed graphs
    if (D::is_directed) {
      NbrsIter vp, vp_end;
      std::tie(vp, vp_end) = boost::adjacent_vertices(vboost, _g);
      for (; vp != vp_end; ++vp) {
        V u = GetV(*vp);
        E e = GetE(boost::edge(*vp, vboost, _g).first);
        _em.erase(e);
        _e.erase(std::find(_e.begin(), _e.end(), e));
      }
    }
    // Remove the vertex
    _vm.erase(v);
    _v.erase(std::find(_v.begin(), _v.end(), v));
    _pre.erase(v);
    if (D::is_directed) _suc.erase(v);
    boost::clear_vertex(vboost, _g);
    boost::remove_vertex(vboost, _g);
  }
  
  BASEGRAPH(void)::AddEdge(const V& u, const V& v, const E& e) {
    ModificationMade();
    if (!HasVertex(u)) AddVertex(u);
    if (!HasVertex(v)) AddVertex(v);
    VertType uboost = GetDescriptor(u);
    VertType vboost = GetDescriptor(v);
    EdgeType eboost = boost::add_edge(uboost, vboost, EP(), _g).first;
    _em.insert(e, eboost);
    _e.emplace_back(e);
  }
  
  BASEGRAPH(void)::RemoveEdge(const E& e) {
    ModificationMade();
    EdgeType eboost = GetDescriptor(e);
    _em.erase(eboost);
    _e.erase(std::find(_e.begin(), _e.end(), e));
    boost::remove_edge(eboost, _g);
  }
  
  BASEGRAPH(void)::RemoveEdge(const V& u, const V& v) {
    E e = GetEdge(u, v);
    RemoveEdge(e);
  }
  
  BASEGRAPH(bool)::HasVertex(const V& v) const {
    return _vm.left.find(v) != _vm.left.end();
  }
  
  BASEGRAPH(bool)::HasEdge(const E& e) const {
    return _em.left.find(e) != _em.left.end();
  }
  
  BASEGRAPH(bool)::HasEdge(const V& u, const V& v) const {
    if (!HasVertex(u) || !HasVertex(v)) return false;
    VertType u_ = GetDescriptor(u);
    VertType v_ = GetDescriptor(v);
    return boost::edge(u_, v_, _g).second;
  }
  
  BASEGRAPH(int64_t)::NumVertices() const {
    return boost::num_vertices(_g);
  }
  
  BASEGRAPH(int64_t)::NumEdges() const {
    return boost::num_edges(_g);
  }
  
  BASEGRAPH(int64_t)::Degree(const V &v) const {
    return HasVertex(v) ? OutDegree(GetDescriptor(v)) : -1;
  }
  
  BASEGRAPH(int64_t)::InDegree(const V &v) const {
    if (!HasVertex(v)) return -1;
    if (D::is_directed) return InDegree(GetDescriptor(v));
    return OutDegree(GetDescriptor(v));
  }
  
  BASEGRAPH(const tBG::VertContain&)::GetNeighbours(const V& v) {
    static_assert(!D::is_directed, "Requires an undirected graph.");
    VertType vboost = GetDescriptor(v);
    auto adjis = boost::adjacent_vertices(vboost, _g);
    _pre.at(v).clear();
    for (; adjis.first != adjis.second; ++adjis.first)
      _pre.at(v).emplace_back(GetV(*adjis.first));
    return _pre.at(v);
  }
  
  BASEGRAPH(const tBG::VertContain&)::GetPredecessors(const V& v) {
    static_assert(D::is_directed, "Requires a directed graph.");
    VertType vboost = GetDescriptor(v);
    auto adjis = boost::inv_adjacent_vertices(vboost, _g);
    _pre.at(v).clear();
    for (; adjis.first != adjis.second; ++adjis.first)
      _pre.at(v).emplace_back(GetV(*adjis.first));
    return _pre.at(v);
  }
  
  BASEGRAPH(const tBG::VertContain&)::GetSuccessors(const V& v) {
    static_assert(D::is_directed, "Requires a directed graph.");
    VertType vboost = GetDescriptor(v);
    auto adjis = boost::adjacent_vertices(vboost, _g);
    _suc.at(v).clear();
    for (; adjis.first != adjis.second; ++adjis.first)
      _suc.at(v).emplace_back(GetV(*adjis.first));
    return _suc.at(v);
  }
  
  BASEGRAPH(std::pair<V,V>)::GetVertices(const E& e) const {
    EdgeType eboost = GetDescriptor(e);
    VertType u = boost::source(eboost, _g);
    VertType v = boost::target(eboost, _g);
    return {GetV(u), GetV(v)};
  }
  
  BASEGRAPH(const tBG::VertContain&)::GetVertices() const {
    return _v;
  }
  
  BASEGRAPH(const tBG::EdgeContain&)::GetEdges() const {
    return _e;
  }
  
  BASEGRAPH(E)::GetEdge(const V& u, const V& v) const {
    VertType uboost = GetDescriptor(u);
    VertType vboost = GetDescriptor(v);
    EdgeType eboost = boost::edge(uboost, vboost, _g).first;
    return GetE(eboost);
  }
  
  BASEGRAPH(V)::GetSourceVertex(const E& e) const {
    EdgeType eboost = GetDescriptor(e);
    VertType source = boost::source(eboost, _g);
    return GetV(source);
  }
  
  BASEGRAPH(V)::GetTargetVertex(const E& e) const {
    EdgeType eboost = GetDescriptor(e);
    VertType target = boost::target(eboost, _g);
    return GetV(target);
  }

// ============================================================================
// == Private helper functions ================================================
// ============================================================================
  BASEGRAPH(tBG::VertType)::GetDescriptor(const V& v) const {
    return _vm.left.at(v);
  }
  
  BASEGRAPH(V)::GetV(BG::VertType v) const {
    return _vm.right.at(v);
  }
  
  BASEGRAPH(tBG::EdgeType)::GetDescriptor(const E& e) const {
    return _em.left.at(e);
  }
  
  BASEGRAPH(E)::GetE(BG::EdgeType e) const {
    return _em.right.at(e);
  }
  
  BASEGRAPH(int64_t)::OutDegree(BG::VertType v) const {
    return boost::out_degree(v, _g);
  }
  
  BASEGRAPH(int64_t)::InDegree(BG::VertType v) const {
    return boost::in_degree(v, _g);
  }
  
// ============================================================================
// == Algorithmic stuffs ======================================================
// ============================================================================
  
  // Connectivity
  BASEGRAPH(bool)::IsConnected() {
    GetConnectedComponents();
    return _comp_cache.size() == 1;
  }
  
  BASEGRAPH(int64_t)::NumConnectedComponents() {
    GetConnectedComponents();
    return _comp_cache.size();
  }

  BASEGRAPH(const tBG::ComponentContain&)::GetConnectedComponents() {
    if (_comp_state == GetCurrentState()) return _comp_cache;
    algorithm::ConnectedComponents(*this, _comp_cache);
    _comp_state = GetCurrentState();
    return _comp_cache;
  }
  
  // Cycles
  BASEGRAPH(bool)::IsCyclic(const V& v) {
    GetCycles();
    return _vcyclic_cache.find(v) != _vcyclic_cache.end();
  }
  
  BASEGRAPH(bool)::IsCyclic(const E& e) {
    GetCycles();
    return _ecyclic_cache.find(e) != _ecyclic_cache.end();
  }
  
  BASEGRAPH(const tBG::CycleEdgeContain&)::GetCycles() {
    using VertSet = eastl::vector_set<V>;
    using EdgeSet = eastl::vector_set<E>;
    
    if (GetCurrentState() == _cycle_state) return _cycles_cache;
    algorithm::AllCycles(*this, _cycles_cache);
    VertSet cyclic_v; cyclic_v.reserve(NumVertices());
    EdgeSet cyclic_e; cyclic_e.reserve(NumEdges());
    
    for (EdgeContain& cycle : _cycles_cache) {
      for (E& edge : cycle) {
        cyclic_e.emplace(edge);
        cyclic_v.emplace(GetSourceVertex(edge));
        cyclic_v.emplace(GetTargetVertex(edge));
      }
    }
    _vcyclic_cache.clear();
    _vcyclic_cache.assign(cyclic_v.begin(), cyclic_v.end());
    _ecyclic_cache.clear();
    _ecyclic_cache.assign(cyclic_e.begin(), cyclic_e.end());
    std::sort(_cycles_cache.begin(), _cycles_cache.end(),
              [](EdgeContain& a, EdgeContain& b){ return a.size() < b.size(); });
    _cycle_state = GetCurrentState();
    return _cycles_cache;
  }
  
  BASEGRAPH(int64_t)::NumCycles() {
    GetCycles();
    return _cycles_cache.size();
  }
  
  
  
  // Isomorphism
  
#undef GRAPHTEMPLATE
#undef BASEGRAPH
#undef BG
#undef tBG
}
