#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <indigox/python/interface.hpp>

#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>

namespace py = pybind11;

void GeneratePyElectronAssignmentGraph(py::module& m) {
  using namespace indigox::graph;
  
  py::class_<IXAGVertex, AGVertex>(m, "AGVertex")
  // Getters
  .def("GetPreAssignedCount", &IXAGVertex::GetPreAssignedCount)
  .def("GetSourceEdge", &IXAGVertex::GetSourceEdge)
  .def("GetSourceVertex", &IXAGVertex::GetSourceVertex)
  .def("GetTotalAssigned", &IXAGVertex::GetTotalAssigned)
  // Checkers
  .def("IsEdgeMapped", &IXAGVertex::IsEdgeMapped)
  .def("IsVertexMapped", &IXAGVertex::IsVertexMapped)
  // Setters
  .def("SetAssignedCount", &IXAGVertex::SetAssignedCount)
  .def("SetPreAssignedCount", &IXAGVertex::SetPreAssignedCount)
  ;
  
  auto get_neighbours = [](const AssignmentGraph& g, const AGVertex& v) {
    auto iters = g->GetNeighbours(v);
    return py::make_iterator(iters.first, iters.second);
  };
  
  auto get_vertices = [](const AssignmentGraph& g) {
    auto iters = g->GetVertices();
    return py::make_iterator(iters.first, iters.second);
  };
  
  py::class_<IXAssignmentGraph, AssignmentGraph>(m, "AssignmentGraph")
  // Constructor
  .def(py::init<MolecularGraph>())
  // Getters
  .def("GetNeighbours", get_neighbours, py::keep_alive<0, 1>())
  .def("GetVertex", py::overload_cast<const MGVertex&>(&IXAssignmentGraph::GetVertex, py::const_))
  .def("GetVertex", py::overload_cast<const MGEdge&>(&IXAssignmentGraph::GetVertex, py::const_))
  .def("GetVertices", get_vertices, py::keep_alive<0, 1>())
  .def("NumVertices", &IXAssignmentGraph::NumVertices)
  // Checkers
  .def("HasVertex", py::overload_cast<const AGVertex&>(&IXAssignmentGraph::HasVertex, py::const_))
  .def("HasVertex", py::overload_cast<const MGVertex&>(&IXAssignmentGraph::HasVertex, py::const_))
  .def("HasVertex", py::overload_cast<const MGEdge&>(&IXAssignmentGraph::HasVertex, py::const_))
//  .def("IsConnected", &IXAssignmentGraph::IsConnected)
  ;
}
