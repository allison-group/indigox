#include <indigox/classes/atom.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox::graph {
  IXAssignmentGraph::IXAssignmentGraph(MolecularGraph g)
  : _source(g), _g() {
    // Add all vertices with edges and the associated edges.
    for (auto it = g->GetEdges(); it.first != it.second; ++it.first) {
      MGEdge edge = *it.first;
      MGVertex s, t;
      std::tie(s,t) = g->GetVertices(edge);
      if (_verts_v.find(s) == _verts_v.end())
        _verts_v.emplace(s, AGVertex(std::make_shared<IXAGVertex>(s)));
      if (_verts_v.find(t) == _verts_v.end())
        _verts_v.emplace(t, AGVertex(std::make_shared<IXAGVertex>(t)));
      _verts_e.emplace(edge, AGVertex(std::make_shared<IXAGVertex>(edge)));
      AGVertex u = _verts_v.at(s), v = _verts_v.at(t);
      AGVertex e = _verts_e.at(edge);
      _g.AddEdge(u.get(), e.get(), nullptr);
      _g.AddEdge(e.get(), v.get(), nullptr);
    }
    // Add unconnected vertices
    for (auto it = g->GetVertices(); it.first != it.second; ++it.first) {
      if (_verts_v.find(*it.first) != _verts_v.end()) continue;
      MGVertex vert = *it.first;
      _verts_v.emplace(vert, AGVertex(std::make_shared<IXAGVertex>(vert)));
      _g.AddVertex(_verts_v.at(vert).get());
    }
    
    _verts.reserve(_verts_v.size() + _verts_e.size());
    // Populate all verts and nbrs
    std::vector<IXAGVertex*> nbrs, verts;
    _g.GetVertices(verts);
    for (auto v : verts) {
      _verts.emplace_back(v->shared_from_this());
      _nbrs.emplace(v->shared_from_this(), std::vector<AGVertex>());
      _nbrs.at(v->shared_from_this()).reserve(_g.Degree(v));
      _g.GetNeighbours(v, nbrs);
      for (auto nbr : nbrs)
        _nbrs.at(v->shared_from_this()).emplace_back(nbr->shared_from_this());
    }
  }
}
