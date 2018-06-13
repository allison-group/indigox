#include <indigox/classes/atom.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox::graph {
  IXAssignmentGraph::IXAssignmentGraph(MolecularGraph g)
  : _source(g), _g() {
    _verts.reserve(g->NumEdges() + g->NumVertices());
    
    // Add all vertices with edges and the associated edges.
    for (auto it = g->GetEdges(); it.first != it.second; ++it.first) {
      MGEdge edge = *it.first;
      MGVertex s, t;
      std::tie(s,t) = g->GetVertices(edge);
      AddEdges(s, t, edge);
    }
    // Add unconnected vertices
    for (auto it = g->GetVertices(); it.first != it.second; ++it.first) {
      if (_verts_v.find(*it.first) != _verts_v.end()) continue;
      AddVertex(*it.first);
    }
    
    // Populate all verts and nbrs
    DetermineAllNeighbours();
  }
  
  AGVertex IXAssignmentGraph::AddVertex(MGVertex v) {
    AGVertex va = std::make_shared<IXAGVertex>(v);
    _verts_v.emplace(v, va);
    _verts.emplace_back(va);
    _g.AddVertex(va.get());
    return va;
  }
  
  void IXAssignmentGraph::AddEdges(MGVertex s, MGVertex t, MGEdge e) {
    AGVertex sv, tv;
    if (_verts_v.find(s) == _verts_v.end()) {
      sv = std::make_shared<IXAGVertex>(s);
      _verts.emplace_back(sv);
      _verts_v.emplace(s, sv);
    }
    else sv = _verts_v.at(s);
    if (_verts_v.find(t) == _verts_v.end()) {
      tv = std::make_shared<IXAGVertex>(t);
      _verts.emplace_back(tv);
      _verts_v.emplace(t, tv);
    }
    else tv = _verts_v.at(t);
    AGVertex ev = std::make_shared<IXAGVertex>(e);
    _verts.emplace_back(ev);
    _verts_e.emplace(e, ev);
    
    _g.AddEdge(sv.get(), ev.get(), nullptr);
    _g.AddEdge(ev.get(), tv.get(), nullptr);
  }
  
  void IXAssignmentGraph::DetermineAllNeighbours() {
    _nbrs.clear();
    std::vector<IXAGVertex*> nbrs;
    for (AGVertex v : _verts) {
      _nbrs.emplace(v, std::vector<AGVertex>());
      _nbrs.at(v).reserve(_g.Degree(v.get()));
      _g.GetNeighbours(v.get(), nbrs);
      for (IXAGVertex* nbr : nbrs)
        _nbrs.at(v).emplace_back(nbr->shared_from_this());
    }
  }
  
  
}
