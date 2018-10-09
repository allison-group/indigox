#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <indigox/python/interface.hpp>

#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

namespace py = pybind11;

void GeneratePyCondensedMolecularGraph(py::module& m) {
  using namespace indigox;
  using namespace indigox::graph;
  
  py::class_<IXCMGVertex, CMGVertex> cmgv(m, "CMGVertex");
  cmgv.def(py::init<const MGVertex&, const CondensedMolecularGraph&>())
  .def("GetSource", &IXCMGVertex::GetSource)
  .def("GetGraph", &IXCMGVertex::GetGraph)
  .def("NumContracted", py::overload_cast<>(&IXCMGVertex::NumContracted, py::const_))
  .def("NumContracted", py::overload_cast<IXCMGVertex::ContractedSymmetry>(&IXCMGVertex::NumContracted, py::const_))
  .def("GetIsomorphismMask", &IXCMGVertex::GetIsomorphismMask)
  .def("IsContractedHere", &IXCMGVertex::IsContractedHere)
  .def("GetContractedVertices", &IXCMGVertex::GetContractedVertices)
  ;
  
  py::enum_<IXCMGVertex::ContractedSymmetry>(cmgv, "ContractedSymmetry")
  .value("Hydrogen", IXCMGVertex::ContractedSymmetry::Hydrogen)
  .value("Fluorine", IXCMGVertex::ContractedSymmetry::Fluorine)
  .value("Chlorine", IXCMGVertex::ContractedSymmetry::Chlorine)
  .value("Bromine", IXCMGVertex::ContractedSymmetry::Bromine)
  .value("Iodine", IXCMGVertex::ContractedSymmetry::Iodine)
  ;
  
  py::class_<IXCMGEdge, CMGEdge>(m, "CMGEdge")
  .def(py::init<const MGEdge&, const CondensedMolecularGraph&>())
  .def("GetSource", &IXCMGEdge::GetSource)
  .def("GetGraph", &IXCMGEdge::GetGraph)
  .def("GetIsomorphismMask", &IXCMGEdge::GetIsomorphismMask)
  ;

  using CMG = IXCondensedMolecularGraph;
  py::class_<CMG, CondensedMolecularGraph>(m, "CondensedMolecularGraph")
  .def(py::init<const MolecularGraph&>())
  //.def("InduceSubgraph"
  //.def("Subgraph"
  .def("GetSource", py::overload_cast<>(&CMG::GetSource, py::const_))
  .def("GetSource", py::overload_cast<const CMGEdge>(&CMG::GetSource, py::const_))
  .def("GetTarget", py::overload_cast<const CMGEdge>(&CMG::GetTarget, py::const_))
  .def("Degree", &CMG::Degree)
  .def("GetEdge", py::overload_cast<const CMGVertex, const CMGVertex>(&CMG::GetEdge, py::const_))
  .def("GetEdge", py::overload_cast<const MGEdge>(&CMG::GetEdge, py::const_))
  .def("GetVertex", &CMG::GetVertex)
  .def("GetEdges", [](CondensedMolecularGraph g) {
    auto ab = g->GetEdges();
    return py::make_iterator(ab.first, ab.second); }, py::keep_alive<0, 1>())
  .def("GetNeighbours", [](CondensedMolecularGraph g, const CMGVertex v) {
    auto ab = g->GetNeighbours(v);
    return py::make_iterator(ab.first, ab.second); }, py::keep_alive<0, 1>())
  .def("GetVertices", py::overload_cast<const CMGEdge>(&CMG::GetVertices, py::const_))
  .def("GetVertices", [](CondensedMolecularGraph g) {
    auto ab = g->GetVertices();
    return py::make_iterator(ab.first, ab.second); }, py::keep_alive<0, 1>())
  .def("HasVertex", py::overload_cast<const MGVertex&>(&CMG::HasVertex, py::const_))
  .def("HasVertex", py::overload_cast<const CMGVertex>(&CMG::HasVertex, py::const_))
  .def("HasCondensedVertex", &CMG::HasCondensedVertex)
  .def("HasEdge", py::overload_cast<const MGEdge&>(&CMG::HasEdge, py::const_))
  .def("HasEdge", py::overload_cast<const CMGEdge>(&CMG::HasEdge, py::const_))
  .def("HasEdge", py::overload_cast<const CMGVertex, const CMGVertex>(&CMG::HasEdge, py::const_))
  .def("NumEdges", &CMG::NumEdges)
  .def("NumVertices", &CMG::NumVertices)
  ;
  
  m.def("CondenseMolecularGraph", &CondenseMolecularGraph);
  
}

