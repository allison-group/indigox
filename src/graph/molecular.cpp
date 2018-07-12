#include <algorithm>
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
    _bn2e.emplace(bnd, e);
    _e.emplace_back(e);
    _n[u].emplace_back(v);
    _n[v].emplace_back(u);
    _g.AddEdge(u.get(), v.get(), e.get());
    return e;
  }
  
  MGVertex IXMolecularGraph::AddVertex(const Atom atm) {
    MGVertex v = MGVertex(new IXMGVertex(atm));
    _at2v.emplace(atm, v);
    _v.emplace_back(v);
    _n.emplace(v, std::vector<MGVertex>());
    _g.AddVertex(v.get());
    return v;
  }
  
  void IXMolecularGraph::Clear() {
    _g.Clear();
    _at2v.clear();
    _bn2e.clear();
    _v.clear();
    _e.clear();
    _n.clear();
    _c.clear();
  }
  
  void IXMolecularGraph::RemoveEdge(const MGEdge e) {
    _g.RemoveEdge(e.get());
    _bn2e.erase(e->GetBond());
    _e.erase(std::find(_e.begin(), _e.end(), e));
    MGVertex u = _at2v[e->GetBond()->GetSourceAtom()];
    MGVertex v = _at2v[e->GetBond()->GetSourceAtom()];
    _n[u].erase(std::find(_n[u].begin(), _n[u].end(), v));
    _n[v].erase(std::find(_n[v].begin(), _n[v].end(), u));
  }
  
  void IXMolecularGraph::RemoveEdge(const MGVertex u, const MGVertex v) {
    RemoveEdge(GetEdge(u, v));
  }
  
  void IXMolecularGraph::RemoveVertex(const MGVertex v) {
    std::vector<MGVertex> nbrs(_n[v].begin(), _n[v].end());
    for (MGVertex n : nbrs) RemoveEdge(v, n);
    _g.RemoveVertex(v.get());
    _at2v.erase(v->GetAtom());
    _v.erase(std::find(_v.begin(), _v.end(), v));
    _n.erase(v);
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
    return _bn2e.at(bnd);
  }
  
  MGVertex IXMolecularGraph::GetVertex(const Atom atm) const {
    if (!HasVertex(atm)) return MGVertex();
    return _at2v.at(atm);
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
  
  
  std::pair<IXMolecularGraph::CompIter, IXMolecularGraph::CompIter>
  IXMolecularGraph::GetConnectedComponents() {
    _c.clear();
    std::vector<std::vector<IXMGVertex*>> ptr_components;
    size_ num = _g.ConnectedComponents(ptr_components);
    _c.reserve(num);
    for (auto component : ptr_components) {
      std::vector<MGVertex> c;
      c.reserve(component.size());
      for (IXMGVertex* v : component)
        c.emplace_back(v->shared_from_this());
      _c.emplace_back(c.begin(), c.end());
    }
    return {_c.cbegin(), _c.cend()};
  }
  
}

