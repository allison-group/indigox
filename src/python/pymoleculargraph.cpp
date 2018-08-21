#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <indigox/python/interface.hpp>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/graph/molecular.hpp>

namespace py = pybind11;

void GeneratePyMolecularGraph(py::module& m) {
  using namespace indigox;
  using namespace indigox::graph;
  
  py::class_<IXMGVertex, MGVertex>(m, "MGVertex")
  .def("GetAtom", &IXMGVertex::GetAtom)
  .def("GetGraph", &IXMGVertex::GetGraph)
  ;
  
  py::class_<IXMGEdge, MGEdge>(m, "MGEdge")
  .def("GetBond", &IXMGEdge::GetBond)
  .def("GetGraph", &IXMGEdge::GetGraph)
  ;
  
//  auto get_connected_components = [] (const MolecularGraph& g) {
//    auto iters = g->GetConnectedComponents();
//    return py::make_iterator(iters.first, iters.second);
//  };
  
  auto get_edges = [](const MolecularGraph& g) {
    auto iters = g->GetEdges();
    return py::make_iterator(iters.first, iters.second);
  };
  
  auto get_neighbours = [](const MolecularGraph& g, const MGVertex& v) {
    auto iters = g->GetNeighbours(v);
    return py::make_iterator(iters.first, iters.second);
  };
  
  auto get_vertices = [](const MolecularGraph& g) {
    auto iters = g->GetVertices();
    return py::make_iterator(iters.first, iters.second);
  };
  
  py::class_<IXMolecularGraph, MolecularGraph>(m, "MolecularGraph")
  // No constructor
  // Getters
//  .def("GetConnectedComponents", get_connected_components, py::keep_alive<0, 1>())
  .def("GetEdge", py::overload_cast<const MGVertex, const MGVertex>(&IXMolecularGraph::GetEdge, py::const_))
  .def("GetEdge", py::overload_cast<const Bond>(&IXMolecularGraph::GetEdge, py::const_))
  .def("GetEdges", get_edges, py::keep_alive<0, 1>())
  .def("GetNeighbours", get_neighbours, py::keep_alive<0, 1>())
  .def("GetSource", &IXMolecularGraph::GetSource)
  .def("GetTarget", &IXMolecularGraph::GetTarget)
  .def("GetVertex", &IXMolecularGraph::GetVertex)
  .def("GetVertices", py::overload_cast<const MGEdge>(&IXMolecularGraph::GetVertices, py::const_))
  .def("GetVertices", get_vertices, py::keep_alive<0, 1>())
//  .def("NumConnectedComponents", &IXMolecularGraph::NumConnectedComponents)
  .def("NumEdges", &IXMolecularGraph::NumEdges)
  .def("NumVertices", &IXMolecularGraph::NumVertices)
  // Checkers
  .def("HasEdge", py::overload_cast<const Bond>(&IXMolecularGraph::HasEdge, py::const_))
  .def("HasEdge", py::overload_cast<const MGEdge>(&IXMolecularGraph::HasEdge, py::const_))
  .def("HasEdge", py::overload_cast<const MGVertex, const MGVertex>(&IXMolecularGraph::HasEdge, py::const_))
  .def("HasVertex", py::overload_cast<const Atom>(&IXMolecularGraph::HasVertex, py::const_))
  .def("HasVertex", py::overload_cast<const MGVertex>(&IXMolecularGraph::HasVertex, py::const_))
//  .def("IsConnected", &IXMolecularGraph::IsConnected)
  // Other
  .def("Degree", &IXMolecularGraph::Degree)
  ;
  
  
}

