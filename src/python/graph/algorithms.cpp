#include <vector>

#include <indigox/python/interface.hpp>
#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/cycles.hpp>
#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/algorithm/cherrypicker.hpp>

#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/parameterised.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(indigox::graph::CondensedMolecularGraph::ComponentContain);
PYBIND11_MAKE_OPAQUE(indigox::graph::CondensedMolecularGraph::CycleEdgeContain);
PYBIND11_MAKE_OPAQUE(indigox::graph::MolecularGraph::ComponentContain);
PYBIND11_MAKE_OPAQUE(indigox::graph::MolecularGraph::CycleEdgeContain);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::sCondensedMolecularGraph>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::sMolecularGraph>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::CMGVertex>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::CMGEdge>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::MGVertex>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::MGEdge>);

void GeneratePyGraphAlgorithms(pybind11::module& m) {
  using namespace indigox;
  using namespace indigox::algorithm;
  using namespace indigox::graph;
  using MG = MolecularGraph;
  using MGS = sMolecularGraph;
  using MGV = MGVertex;
  using MGE = MGEdge;
  using VMG = std::vector<sMolecularGraph>;
  using VMGV = MolecularGraph::VertContain;
  using VMGE = MolecularGraph::EdgeContain;
  using VVMGV = MolecularGraph::ComponentContain;
  using VVMGE = MolecularGraph::CycleEdgeContain;
  
  using CMG = CondensedMolecularGraph;
  using CMGS = sCondensedMolecularGraph;
  using CMGV = CMGVertex;
  using CMGE = CMGEdge;
  using VCMG = std::vector<sCondensedMolecularGraph>;
  using VCMGV = CondensedMolecularGraph::VertContain;
  using VCMGE = CondensedMolecularGraph::EdgeContain;
  using VVCMGV = CondensedMolecularGraph::ComponentContain;
  using VVCMGE = CondensedMolecularGraph::CycleEdgeContain;
  
  using D = Directed;
  using GL = GraphLabel;
  
  m.def("ShortestPath", [](MG& g, MGV u, MGV v) { return ShortestPath(g, u, v); });
  m.def("ShortestPath", [](CMG& g, CMGV u, CMGV v) { return ShortestPath(g, u, v); });
  m.def("AllSimplePaths", [](MG& g, MGV u, MGV v, VVMGE& p) { AllSimplePaths(g, u, v, p); });
  m.def("AllSimplePaths", [](MG& g, MGV u, MGV v, VVMGE& p, int64_t sz) { AllSimplePaths(g, u, v, p, sz); });
  m.def("AllSimplePaths", [](CMG& g, CMGV u, CMGV v, VVCMGE& p) { AllSimplePaths(g, u, v, p); });
  m.def("AllSimplePaths", [](CMG& g, CMGV u, CMGV v, VVCMGE& p, int64_t sz) { AllSimplePaths(g, u, v, p, sz); });
  
  m.def("ConnectedComponents", [](MG& g, VVMGV& c) { return ConnectedComponents(g, c); });
  m.def("ConnectedComponents", [](CMG& g, VVCMGV& c) { return ConnectedComponents(g, c); });
  m.def("ConnectedSubgraphs", [](MG& g, VMG& s) { return ConnectedSubgraphs(g, s); });
  m.def("ConnectedSubgraphs", [](MG& g, VMG& s, int64_t m) { return ConnectedSubgraphs(g, s, m); });
  m.def("ConnectedSubgraphs", [](MG& g, VMG& s, int64_t m, int64_t M) { return ConnectedSubgraphs(g, s, m, M); });
  m.def("ConnectedSubgraphs", [](CMG& g, VCMG& s) { return ConnectedSubgraphs(g, s); });
  m.def("ConnectedSubgraphs", [](CMG& g, VCMG& s, int64_t m) { return ConnectedSubgraphs(g, s, m); });
  m.def("ConnectedSubgraphs", [](CMG& g, VCMG& s, int64_t m, int64_t M) { return ConnectedSubgraphs(g, s, m, M); });
  
  m.def("CycleBasis", [](MG& g, VVMGV& b) { return CycleBasis(g, b); });
  m.def("CycleBasis", [](MG& g, VVMGE& b) { return CycleBasis(g, b); });
  m.def("CycleBasis", [](CMG& g, VVCMGV& b) { return CycleBasis(g, b); });
  m.def("CycleBasis", [](CMG& g, VVCMGE& b) { return CycleBasis(g, b); });
  m.def("AllCycles", [](MG& g, VVMGE& c) { return AllCycles(g, c); });
  m.def("AllCycles", [](CMG& g, VVCMGE& c) { return AllCycles(g, c); });
  
  enum class IsomorphismCallbackTypes {
    Print
  };
  
  py::enum_<IsomorphismCallbackTypes>(m, "IsomorphismCallbackType")
  .value("Print", IsomorphismCallbackTypes::Print);
  
  m.def("SubgraphIsomorphisms", [](CMG& small, CMG& large, IsomorphismCallbackTypes type){
    if (type == IsomorphismCallbackTypes::Print) {
      CMGPrintCallback cb;
      SubgraphIsomorphisms(small, large, cb);
    }
  });
//  m.def("SubgraphIsomorphisms", py::overload_cast<MG&, MG&, MGCallback&>(&SubgraphIsomorphisms));
  
  using VParam = CherryPicker::VertexParameters;
  using EParam = CherryPicker::EdgeParameters;
  using CPSet = CherryPicker::Settings;
  
  py::class_<CherryPicker> cp(m, "CherryPicker");
  
  py::enum_<VParam>(cp, "VertexParameters", py::arithmetic())
  .value("None", VParam::None)
  .value("ElementType", VParam::ElementType)
  .value("FormalCharge", VParam::FormalCharge)
  .value("CondensedVertices", VParam::CondensedVertices)
  .value("CyclicNature", VParam::CyclicNature)
  .value("Stereochemistry", VParam::Stereochemistry)
  .value("Aromaticity", VParam::Aromaticity)
  ;
  
  py::enum_<EParam>(cp, "EdgeParameters", py::arithmetic())
  .value("None", EParam::None)
  .value("BondOrder", EParam::BondOrder)
  .value("Stereochemistry", EParam::Stereochemistry)
  .value("CyclicNature", EParam::CyclicNature)
  .value("Aromaticity", EParam::Aromaticity)
  ;
  
  py::class_<CPSet>(cp, "Settings")
  .def_readwrite_static("AllowDanglingBonds", &CPSet::AllowDanglingBonds)
  .def_readwrite_static("AllowDanglingAngles", &CPSet::AllowDanglingAngles)
  .def_readwrite_static("AllowDanglingDihedrals", &CPSet::AllowDanglingDihedrals)
  .def_readwrite_static("VertexMapping", &CPSet::VertexMapping)
  .def_readwrite_static("EdgeMapping", &CPSet::EdgeMapping)
  ;
  
  cp.def(py::init<const Forcefield&>())
  .def("AddAthenaeum", &CherryPicker::AddAthenaeum)
  .def("RemoveAthenaeum", &CherryPicker::RemoveAthenaeum)
  .def("NumAthenaeums", &CherryPicker::NumAthenaeums)
  .def("ParameteriseMolecule", &CherryPicker::ParameteriseMolecule)
  .def("GetForcefield", &CherryPicker::GetForcefield)
  ;
  
}
