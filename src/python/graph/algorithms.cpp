#include <indigox/algorithm/cherrypicker.hpp>
#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/cycles.hpp>
#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/python/interface.hpp>

#include <vector>

namespace py = pybind11;

void GeneratePyGraphAlgorithms(pybind11::module &m) {
  using namespace indigox;
  using namespace indigox::algorithm;
  using namespace indigox::graph;
  using MG = MolecularGraph;
  //  using MGS = MolecularGraph;
  using MGV = MGVertex;
  //  using MGE = MGEdge;
  using VMG = std::vector<MolecularGraph>;
  //  using VMGV = MolecularGraph::VertContain;
  //  using VMGE = MolecularGraph::EdgeContain;
  using VVMGV = MolecularGraph::ComponentContain;
  using VVMGE = MolecularGraph::CycleEdgeContain;

  using CMG = CondensedMolecularGraph;
  //  using CMGS = CondensedMolecularGraph;
  using CMGV = CMGVertex;
  //  using CMGE = CMGEdge;
  using VCMG = std::vector<CondensedMolecularGraph>;
  //  using VCMGV = CondensedMolecularGraph::VertContain;
  //  using VCMGE = CondensedMolecularGraph::EdgeContain;
  using VVCMGV = CondensedMolecularGraph::ComponentContain;
  using VVCMGE = CondensedMolecularGraph::CycleEdgeContain;

  //  using D = Directed;
  //  using GL = GraphLabel;

  m.def("ShortestPath",
        [](MG &g, MGV u, MGV v) { return ShortestPath(g, u, v); },
        py::arg("graph"), py::arg("source"), py::arg("target"));
  m.def("ShortestPath",
        [](CMG &g, CMGV u, CMGV v) { return ShortestPath(g, u, v); },
        py::arg("graph"), py::arg("source"), py::arg("target"));

  m.def("AllSimplePaths",
        [](MG &g, MGV u, MGV v, int64_t sz) -> VVMGE {
          VVMGE p;
          AllSimplePaths(g, u, v, p, sz);
          return p;
        },
        py::arg("graph"), py::arg("source"), py::arg("target"),
        py::arg("maximum_length") = -1);
  m.def("AllSimplePaths",
        [](CMG &g, CMGV u, CMGV v, int64_t sz) -> VVCMGE {
          VVCMGE p;
          AllSimplePaths(g, u, v, p, sz);
          return p;
        },
        py::arg("graph"), py::arg("source"), py::arg("target"),
        py::arg("maximum_length") = -1);

  m.def("ConnectedComponents",
        [](MG &g, VVMGV &c) { return ConnectedComponents(g, c); });
  m.def("ConnectedComponents",
        [](CMG &g, VVCMGV &c) { return ConnectedComponents(g, c); });

  m.def("ConnectedSubgraphs",
        [](MG &g, size_t min, size_t max) -> VMG {
          ConnectedSubgraphs gen(g, min, max);
          MG subg;
          VMG subgraphs;
          while (gen(subg)) { subgraphs.emplace_back(subg); }
          return subgraphs;
        },
        py::arg("graph"), py::arg("minimum_size") = 0,
        py::arg("maximum_size") = std::numeric_limits<size_t>::max());
  m.def("ConnectedSubgraphs",
        [](CMG &g, size_t min, size_t max) -> VCMG {
          ConnectedSubgraphs gen(g, min, max);
          CMG subg;
          VCMG subgraphs;
          while (gen(subg)) { subgraphs.emplace_back(subg); }
          return subgraphs;
        },
        py::arg("graph"), py::arg("minimum_size") = 0,
        py::arg("maximum_size") = std::numeric_limits<size_t>::max());

  m.def("OptimalChargeGroups", &OptimalChargeGroups, py::arg("molecule"),
        py::arg("size_limit") = 5);

  m.def("CycleBasis", [](MG &g, VVMGV &b) { return CycleBasis(g, b); });
  m.def("CycleBasis", [](MG &g, VVMGE &b) { return CycleBasis(g, b); });
  m.def("CycleBasis", [](CMG &g, VVCMGV &b) { return CycleBasis(g, b); });
  m.def("CycleBasis", [](CMG &g, VVCMGE &b) { return CycleBasis(g, b); });
  m.def("AllCycles", [](MG &g, VVMGE &c) { return AllCycles(g, c); });
  m.def("AllCycles", [](CMG &g, VVCMGE &c) { return AllCycles(g, c); });

  enum class IsomorphismCallbackTypes { Print };

  py::enum_<IsomorphismCallbackTypes>(m, "IsomorphismCallbackType")
      .value("Print", IsomorphismCallbackTypes::Print);

  m.def("SubgraphIsomorphisms",
        [](CMG &small, CMG &large, IsomorphismCallbackTypes type) {
          if (type == IsomorphismCallbackTypes::Print) {
            CMGPrintCallback cb;
            SubgraphIsomorphisms(small, large, cb);
          }
        });
  //  m.def("SubgraphIsomorphisms", py::overload_cast<MG&, MG&,
  //  MGCallback&>(&SubgraphIsomorphisms));

  using CPSet = CherryPicker::Settings;

  py::class_<CherryPicker> cherrypicker(m, "CherryPicker");

  py::enum_<CPSet>(cherrypicker, "Settings")
      // Boolean settings
      .value("VertexElement", CPSet::VertexElement)
      .value("VertexFormalCharge", CPSet::VertexFormalCharge)
      .value("VertexCondensed", CPSet::VertexCondensed)
      .value("VertexCyclic", CPSet::VertexCyclic)
      .value("VertexCyclicSize", CPSet::VertexCyclicSize)
      .value("VertexAromaticity", CPSet::VertexAromaticity)
      .value("VertexDegree", CPSet::VertexDegree)
      .value("EdgeBondOrder", CPSet::EdgeBondOrder)
      .value("EdgeStereochemistry", CPSet::EdgeStereochemistry)
      .value("EdgeCyclic", CPSet::EdgeCyclic)
      .value("EdgeCyclicSize", CPSet::EdgeCyclicSize)
      .value("EdgeDegree", CPSet::EdgeDegree)
      .value("AllowDanglingBonds", CPSet::AllowDanglingBonds)
      .value("AllowDanglingAngles", CPSet::AllowDanglingAngles)
      .value("AllowDanglingDihedrals", CPSet::AllowDanglingDihedrals)
      .value("ParameteriseFromAllPermutations",
             CPSet::ParameteriseFromAllPermutations)
      // Integer settings
      .value("MinimumFragmentSize", CPSet::MinimumFragmentSize)
      .value("MaximumFragmentSize", CPSet::MaximumFragmentSize);

  cherrypicker.def(py::init<Forcefield &>())
      .def("AddAthenaeum", &CherryPicker::AddAthenaeum)
      .def("RemoveAthenaeum", &CherryPicker::RemoveAthenaeum)
      .def("NumAthenaeums", &CherryPicker::NumAthenaeums)
      .def("ParameteriseMolecule", &CherryPicker::ParameteriseMolecule)
      .def("GetForcefield", &CherryPicker::GetForcefield)
      .def("GetBool", &CherryPicker::GetBool)
      .def("SetBool", &CherryPicker::SetBool)
      .def("UnsetBool", &CherryPicker::UnsetBool)
      .def("GetInt", &CherryPicker::GetInt)
      .def("SetInt", &CherryPicker::SetInt)
      .def("DefaultSettings", &CherryPicker::DefaultSettings);
}
