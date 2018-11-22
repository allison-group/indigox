#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

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

void GeneratePyGraphs(py::module& m) {
  using namespace indigox;
  using namespace indigox::graph;
  py::return_value_policy Ref = py::return_value_policy::reference;
  // ===========================================================================
  // == Enum Types bindings ====================================================
  // ===========================================================================
  py::enum_<ContractedSymmetry>(m, "ContractedSymmetry")
  .value("Hydrogen", ContractedSymmetry::Hydrogen)
  .value("Fluorine", ContractedSymmetry::Fluorine)
  .value("Chlorine", ContractedSymmetry::Chlorine)
  .value("Bromine", ContractedSymmetry::Bromine)
  .value("Iodine", ContractedSymmetry::Iodine)
  ;
  
  // ===========================================================================
  // == MGVertex class bindings ================================================
  // ===========================================================================
  py::class_<MGVertex>(m, "MGVertex")
  .def("GetAtom", &MGVertex::GetAtom, Ref)
  .def("GetGraph", &MGVertex::GetGraph, Ref)
  .def(py::init<>())
  .def(py::init<const MGVertex&>())
  .def(py::self == py::self)
  .def(py::self < py::self)
  ;
  
  // ===========================================================================
  // == MGEdge class bindings ==================================================
  // ===========================================================================
  py::class_<MGEdge>(m, "MGEdge")
  .def("GetBond", &MGEdge::GetBond, Ref)
  .def("GetGraph", &MGEdge::GetGraph, Ref)
  .def(py::init<>())
  .def(py::init<const MGEdge&>())
  .def(py::self == py::self)
  .def(py::self < py::self)
  ;
  
  // ===========================================================================
  // == MolecularGraph class bindings ==========================================
  // ===========================================================================
  using MGV = MGVertex;
  using MGE = MGEdge;
  using MG = MolecularGraph;
  py::class_<MG, sMolecularGraph>(m, "MolecularGraph")
  .def("Subgraph", py::overload_cast<std::vector<MGV>&>(&MG::Subgraph))
  .def("Subgraph", py::overload_cast<std::vector<MGV>&, std::vector<MGE>&>(&MG::Subgraph))
  .def("IsSubgraph", &MG::IsSubgraph)
  .def("GetEdge", py::overload_cast<Bond&>(&MG::GetEdge, py::const_))
  .def("GetVertex", &MG::GetVertex)
  .def("HasVertex", py::overload_cast<Atom&>(&MG::HasVertex, py::const_))
  .def("HasEdge", py::overload_cast<Bond&>(&MG::HasEdge, py::const_))
  .def("GetMolecule", &MG::GetMolecule, Ref)
  .def("GetSuperGraph", &MG::GetSuperGraph, Ref)
  .def("GetCondensedGraph", &MG::GetCondensedGraph, Ref)
  
  .def("HasVertex", py::overload_cast<const MGV&>(&MG::HasVertex, py::const_))
  .def("HasEdge", py::overload_cast<const MGEdge&>(&MG::HasEdge, py::const_))
  .def("HasEdge", py::overload_cast<const MGV&, const MGV&>(&MG::HasEdge, py::const_))
  .def("NumVertices", &MG::NumVertices)
  .def("NumEdges", &MG::NumEdges)
  .def("Degree", &MG::Degree)
  .def("InDegree", py::overload_cast<const MGV&>(&MG::InDegree, py::const_))
  .def("GetNeighbours", &MG::GetNeighbours, Ref)
  .def("GetVertices", py::overload_cast<const MGE&>(&MG::GetVertices, py::const_))
  .def("GetVertices", py::overload_cast<>(&MG::GetVertices, py::const_))
  .def("GetEdges", py::overload_cast<>(&MG::GetEdges, py::const_))
  .def("GetEdge", py::overload_cast<const MGV&, const MGV&>(&MG::GetEdge, py::const_))
  .def("GetSourceVertex", &MG::GetSourceVertex)
  .def("GetTargetVertex", &MG::GetTargetVertex)
  .def("IsConnected", &MG::IsConnected)
  .def("NumConnectedComponents", &MG::NumConnectedComponents)
  .def("GetConnectedComponents", &MG::GetConnectedComponents, Ref)
  .def("IsCyclic", py::overload_cast<const MGV&>(&MG::IsCyclic))
  .def("IsCyclic", py::overload_cast<const MGV&, uint32_t>(&MG::IsCyclic))
  .def("IsCyclic", py::overload_cast<const MGE&>(&MG::IsCyclic))
  .def("IsCyclic", py::overload_cast<const MGE&, uint32_t>(&MG::IsCyclic))
  .def("GetCycles", &MG::GetCycles, Ref)
  .def("NumCycles", &MG::NumCycles)
  ;
  
  // ===========================================================================
  // == CMGVertex class bindings ===============================================
  // ===========================================================================
  py::class_<CMGVertex>(m, "CMGVertex")
  .def(py::init<>())
  .def(py::init<const CMGVertex&>())
  .def(py::self == py::self)
  .def(py::self < py::self)
  .def("GetSource", &CMGVertex::GetSource)
  .def("GetGraph", &CMGVertex::GetGraph, Ref)
  .def("NumContracted", py::overload_cast<>(&CMGVertex::NumContracted, py::const_))
  .def("NumContracted", py::overload_cast<ContractedSymmetry>(&CMGVertex::NumContracted, py::const_))
  .def("GetIsomorphismMask", &CMGVertex::GetIsomorphismMask)
  .def("IsContractedHere", &CMGVertex::IsContractedHere)
  .def("GetContractedVertices", &CMGVertex::GetContractedVertices)
  .def("GetCondensedVertices", &CMGVertex::GetCondensedVertices)
  ;
  
  // ===========================================================================
  // == CMGEdge class bindings =================================================
  // ===========================================================================
  py::class_<CMGEdge>(m, "CMGEdge")
  .def(py::init<>())
  .def(py::init<const CMGEdge&>())
  .def(py::self == py::self)
  .def(py::self < py::self)
  .def("GetSource", &CMGEdge::GetSource)
  .def("GetGraph", &CMGEdge::GetGraph)
  .def("GetIsomorphismMask", &CMGEdge::GetIsomorphismMask)
  ;
  
  // ===========================================================================
  // == CondensedMolecularGraph class bindings =================================
  // ===========================================================================
  using CMGV = CMGVertex;
  using CMGE = CMGEdge;
  using CMG = CondensedMolecularGraph;
  py::class_<CMG, sCondensedMolecularGraph>(m, "CondensedMolecularGraph")
  .def("Subgraph", py::overload_cast<std::vector<CMGV>&>(&CMG::Subgraph))
  .def("Subgraph", py::overload_cast<std::vector<CMGV>&, std::vector<CMGE>&>(&CMG::Subgraph))
  .def("IsSubgraph", &CMG::IsSubgraph)
  .def("GetMolecularGraph", &CMG::GetMolecularGraph, Ref)
  .def("GetSuperGraph", &CMG::GetSuperGraph, Ref)
  .def("GetEdge", py::overload_cast<const MGE&>(&CMG::GetEdge, py::const_))
  .def("GetVertex", &CMG::GetVertex)
  .def("GetCondensedVertex", &CMG::GetCondensedVertex)
  .def("HasVertex", py::overload_cast<const MGV&>(&CMG::HasVertex, py::const_))
  .def("HasCondensedVertex", &CMG::HasCondensedVertex)
  .def("HasEdge", py::overload_cast<const MGE&>(&CMG::HasEdge, py::const_))
  
  .def("HasVertex", py::overload_cast<const CMGV&>(&CMG::HasVertex, py::const_))
  .def("HasEdge", py::overload_cast<const CMGEdge&>(&CMG::HasEdge, py::const_))
  .def("HasEdge", py::overload_cast<const CMGV&, const CMGV&>(&CMG::HasEdge, py::const_))
  .def("NumVertices", &CMG::NumVertices)
  .def("NumEdges", &CMG::NumEdges)
  .def("Degree", &CMG::Degree)
  .def("InDegree", py::overload_cast<const CMGV&>(&CMG::InDegree, py::const_))
  .def("GetNeighbours", &CMG::GetNeighbours, Ref)
  .def("GetVertices", py::overload_cast<const CMGE&>(&CMG::GetVertices, py::const_))
  .def("GetVertices", py::overload_cast<>(&CMG::GetVertices, py::const_))
  .def("GetEdges", py::overload_cast<>(&CMG::GetEdges, py::const_))
  .def("GetEdge", py::overload_cast<const CMGV&, const CMGV&>(&CMG::GetEdge, py::const_))
  .def("GetSourceVertex", &CMG::GetSourceVertex)
  .def("GetTargetVertex", &CMG::GetTargetVertex)
  .def("IsConnected", &CMG::IsConnected)
  .def("NumConnectedComponents", &CMG::NumConnectedComponents)
  .def("GetConnectedComponents", &CMG::GetConnectedComponents, Ref)
  .def("IsCyclic", py::overload_cast<const CMGV&>(&CMG::IsCyclic))
  .def("IsCyclic", py::overload_cast<const CMGV&, uint32_t>(&CMG::IsCyclic))
  .def("IsCyclic", py::overload_cast<const CMGE&>(&CMG::IsCyclic))
  .def("IsCyclic", py::overload_cast<const CMGE&, uint32_t>(&CMG::IsCyclic))
  .def("GetCycles", &CMG::GetCycles, Ref)
  .def("NumCycles", &CMG::NumCycles)
  ;
  
  py::bind_vector<CondensedMolecularGraph::VertContain>(m, "VectorCMGVertex");
  py::bind_vector<CondensedMolecularGraph::EdgeContain>(m, "VectorCMGEdge");
  py::bind_vector<MolecularGraph::VertContain>(m, "VectorMGVertex");
  py::bind_vector<MolecularGraph::EdgeContain>(m, "VectorMGEdge");
  py::bind_vector<CondensedMolecularGraph::ComponentContain>(m, "VectorVectorCMGVertex");
  py::bind_vector<CondensedMolecularGraph::CycleEdgeContain>(m, "VectorVectorCMGEdge");
  py::bind_vector<MolecularGraph::ComponentContain>(m, "VectorVectorMGVertex");
  py::bind_vector<MolecularGraph::CycleEdgeContain>(m, "VectorVectorMGEdge");
  py::bind_vector<std::vector<sCondensedMolecularGraph>>(m, "VectorCondensedMolecularGraph");
  py::bind_vector<std::vector<sMolecularGraph>>(m, "VectorMolecularGraph");
  
}

