#include <indigox/classes/atom.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>

namespace indigox::graph {
  IXAssignmentGraph::IXAssignmentGraph(MolecularGraph g)
  : _source(g), _g() {
    _v.reserve(g->NumEdges() + g->NumVertices());
    
    // Add all vertices with edges and the associated edges.
    for (auto it = g->GetEdges(); it.first != it.second; ++it.first) {
      MGEdge edge = *it.first;
      MGVertex s, t;
      std::tie(s,t) = g->GetVertices(edge);
      AddEdges(s, t, edge);
    }
    // Add unconnected vertices
    for (auto it = g->GetVertices(); it.first != it.second; ++it.first) {
      if (_v2v.find(*it.first) != _v2v.end()) continue;
      AddVertex(*it.first);
    }
    
    // Populate all verts and nbrs
    DetermineAllNeighbours();
  }
  
  AGVertex IXAssignmentGraph::AddVertex(MGVertex v) {
    AGVertex va = AGVertex(new IXAGVertex(v));
    _v2v.emplace(v, va);
    _v.emplace_back(va);
    _g.AddVertex(va.get());
    return va;
  }
  
  void IXAssignmentGraph::AddEdges(MGVertex s, MGVertex t, MGEdge e) {
    AGVertex sv, tv;
    if (_v2v.find(s) == _v2v.end()) {
      sv = AGVertex(new IXAGVertex(s));
      _v.emplace_back(sv);
      _v2v.emplace(s, sv);
    }
    else sv = _v2v.at(s);
    if (_v2v.find(t) == _v2v.end()) {
      tv = AGVertex(new IXAGVertex(t));
      _v.emplace_back(tv);
      _v2v.emplace(t, tv);
    }
    else tv = _v2v.at(t);
    AGVertex ev = AGVertex(new IXAGVertex(e));
    ev->SetPreAssignedCount(2);
    _v.emplace_back(ev);
    _e2v.emplace(e, ev);
    
    _g.AddEdge(sv.get(), ev.get(), nullptr);
    _g.AddEdge(ev.get(), tv.get(), nullptr);
  }
  
  void IXAssignmentGraph::DetermineAllNeighbours() {
    _n.clear();
    std::vector<IXAGVertex*> nbrs;
    for (AGVertex v : _v) {
      _n.emplace(v, std::vector<AGVertex>());
      _n.at(v).reserve(_g.Degree(v.get()));
      _g.GetNeighbours(v.get(), nbrs);
      for (IXAGVertex* nbr : nbrs)
        _n.at(v).emplace_back(nbr->shared_from_this());
    }
  }
  
  uint32_t __VertexPrePlace(size_t degree, const Element& element) {
    switch (degree) {
      case 1:
        switch (element.GetAtomicNumber()) {
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
        switch (element.GetAtomicNumber()) {
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
    for (AGVertex v : _v) {
      uint32_t preassign = 0;
      if (v->IsEdgeMapped()) {
        preassign = 2; // no preassigning is performed on bonds
      } else {
        preassign = __VertexPrePlace(Degree(v),
                                v->GetSourceVertex()->GetAtom()->GetElement());
      }
      v->SetPreAssignedCount(preassign);
    }
  }
  
}
