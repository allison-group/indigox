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
    AGVertex va = AGVertex(new IXAGVertex(v));
    _verts_v.emplace(v, va);
    _verts.emplace_back(va);
    _g.AddVertex(va.get());
    return va;
  }
  
  void IXAssignmentGraph::AddEdges(MGVertex s, MGVertex t, MGEdge e) {
    AGVertex sv, tv;
    if (_verts_v.find(s) == _verts_v.end()) {
      sv = AGVertex(new IXAGVertex(s));
      _verts.emplace_back(sv);
      _verts_v.emplace(s, sv);
    }
    else sv = _verts_v.at(s);
    if (_verts_v.find(t) == _verts_v.end()) {
      tv = AGVertex(new IXAGVertex(t));
      _verts.emplace_back(tv);
      _verts_v.emplace(t, tv);
    }
    else tv = _verts_v.at(t);
    AGVertex ev = AGVertex(new IXAGVertex(e));
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
  
  uint_ __VertexPrePlace(size_ degree, const Element& element) {
    switch (degree) {
      case 1:
        switch (element->GetAtomicNumber()) {
          case 9:   // F
          case 17:  // Cl
          case 35:  // Br
            return 6;
          case 8:   // O
          case 16:  // S
            return 4;
          case 7:   // N
            return 2;
          default:  // All others
            return 0;
        }
        break;
      case 2:
        switch (element->GetAtomicNumber()) {
          case 8:   // O
          case 16:  // S
            return 4;
          default:
            return 0;
        }
      default:
        return 0;
    }
  }
  
  void IXAssignmentGraph::PreassignElectrons() {
    for (AGVertex v : _verts) {
      uint_ preassign = 0;
      if (v->IsEdgeMapped()) {  // no preassigning is performed on bonds
      } else {
        preassign = __VertexPrePlace(Degree(v),
                                v->GetSourceVertex()->GetAtom()->GetElement());
      }
      v->SetPreAssignedCount(preassign);
    }
  }
  
}
