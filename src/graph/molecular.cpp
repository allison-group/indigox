#include <limits>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox::graph {
  
  MGEdge IXMolecularGraph::AddEdge(const Bond bnd) {
    Atom source = bnd->GetSourceAtom();
    Atom target = bnd->GetTargetAtom();
    if (!HasVertex(source)) AddVertex(source);
    if (!HasVertex(target)) AddVertex(target);
    // Edges will only be added when a bond is added to the molecule
    MGVertex u = GetVertex(source);
    MGVertex v = GetVertex(target);
    MGEdge e = MGEdge(new IXMGEdge(bnd));
    _edges.emplace(bnd, e);
    _g.AddEdge(u.get(), v.get(), e.get());
    return e;
  }
  
  MGVertex IXMolecularGraph::AddVertex(const Atom atm) {
    _verts.emplace(atm, MGVertex(new IXMGVertex(atm)));
    _g.AddVertex(_verts.at(atm).get());
    return _verts.at(atm);
  }
  
  void IXMolecularGraph::Clear() {
    _g.Clear();
    _verts.clear();
    _edges.clear();
  }
  
  void IXMolecularGraph::RemoveEdge(const MGEdge e) {
    _g.RemoveEdge(e.get());
    _edges.erase(e->GetBond());
  }
  
  void IXMolecularGraph::RemoveEdge(const MGVertex u, const MGVertex v) {
    RemoveEdge(GetEdge(u, v));
  }
  
  void IXMolecularGraph::RemoveVertex(const MGVertex v) {
    std::vector<IXMGVertex*> nbrs;
    _g.GetNeighbours(v.get(), nbrs);
    for (auto n : nbrs) RemoveEdge(v, n->shared_from_this());
    _g.RemoveVertex(v.get());
    _verts.erase(v->GetAtom());
  }
  
  size_ IXMolecularGraph::Degree(const MGVertex v) const {
    if (!HasVertex(v)) return std::numeric_limits<size_>::max();
    return _g.Degree(v.get());
  }
  
  MGEdge IXMolecularGraph::GetEdge(const MGVertex u, const MGVertex v) const {
    if (!HasEdge(u, v)) return MGEdge();
    return _g.GetEdge(u.get(), v.get())->shared_from_this();
  }
  
  MGEdge IXMolecularGraph::GetEdge(const Bond bnd) const {
    if (!HasEdge(bnd)) return MGEdge();
    return _edges.at(bnd);
  }
  
  MGVertex IXMolecularGraph::GetVertex(const Atom atm) const {
    if (!HasVertex(atm)) return MGVertex();
    return _verts.at(atm);
  }
  
  MGVertex IXMolecularGraph::GetSource(const MGEdge e) const {
    if (!HasEdge(e)) return MGVertex();
    return _g.GetSource(e.get())->shared_from_this();
  }
  
  MGVertex IXMolecularGraph::GetTarget(const MGEdge e) const {
    if (!HasEdge(e)) return MGVertex();
    return _g.GetTarget(e.get())->shared_from_this();
  }
  
  std::pair<MGVertex, MGVertex>
  IXMolecularGraph::GetVertices(const MGEdge e) const {
    if (!HasEdge(e)) return {MGVertex(), MGVertex()};
    auto res = _g.GetVertices(e.get());
    return {res.first->shared_from_this(), res.second->shared_from_this()};
  }
  
  std::pair<IXMolecularGraph::VertIter, IXMolecularGraph::VertIter>
  IXMolecularGraph::GetVertices() {
    _vert_access.clear(); _vert_access.reserve(NumVertices());
    std::vector<IXMGVertex*> vert_itrs; _g.GetVertices(vert_itrs);
    for (auto& i : vert_itrs) _vert_access.emplace_back(i->shared_from_this());
    return {_vert_access.cbegin(), _vert_access.cend()};
  }
  
  std::pair<IXMolecularGraph::EdgeIter, IXMolecularGraph::EdgeIter>
  IXMolecularGraph::GetEdges() {
    _edge_access.clear(); _edge_access.reserve(NumEdges());
    std::vector<IXMGEdge*> edge_itrs; _g.GetEdges(edge_itrs);
    for (auto i : edge_itrs) _edge_access.emplace_back(i->shared_from_this());
    return {_edge_access.cbegin(), _edge_access.cend()};
  }
  
  std::pair<IXMolecularGraph::NbrsIter, IXMolecularGraph::NbrsIter>
  IXMolecularGraph::GetNeighbours(const MGVertex v) {
    _nbrs_access.clear();
    if (HasVertex(v)) {
      std::vector<IXMGVertex*> nbrs; _g.GetNeighbours(v.get(), nbrs);
      _nbrs_access.reserve(nbrs.size());
      for (auto i : nbrs) _nbrs_access.emplace_back(i->shared_from_this());
    }
    return {_nbrs_access.cbegin(), _nbrs_access.cend()};
  }
  
  std::pair<IXMolecularGraph::CompIter, IXMolecularGraph::CompIter>
  IXMolecularGraph::GetConnectedComponents() {
    _components.clear();
    std::vector<std::vector<IXMGVertex*>> ptr_components;
    size_ num = _g.ConnectedComponents(ptr_components);
    _components.reserve(num);
    for (auto component : ptr_components) {
      _components.push_back(std::vector<MGVertex>());
      _components.back().reserve(component.size());
      for (IXMGVertex* v : component)
        _components.back().emplace_back(v->shared_from_this());
    }
    return {_components.cbegin(), _components.cend()};
  }
  
}
