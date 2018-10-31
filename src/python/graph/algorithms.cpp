#include <vector>

#include <pybind11/pybind11.h>

#include <indigox/python/interface.hpp>
#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/cycles.hpp>
#include <indigox/algorithm/graph/paths.hpp>

#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

namespace py = pybind11;

void GeneratePyGraphAlgorithms(pybind11::module& m) {
  using namespace indigox::algorithm;
  using namespace indigox::graph;
  using V1 = CMGVertex;
  using E1 = CMGEdge;
  using S1 = sCondensedMolecularGraph;
  using V2 = MGVertex;
  using E2 = MGEdge;
  using S2 = sMolecularGraph;
  using D = Undirected;
  using P = GraphLabel;
  
  using CCContain1 = CondensedMolecularGraph::ComponentContain;
  using CCContain2 = MolecularGraph::ComponentContain;
  m.def("ConnectedComponents", &ConnectedComponents<V1,E1,S1,D,P,P,CCContain1>);
  m.def("ConnectedComponents", &ConnectedComponents<V2,E2,S2,D,P,P,CCContain2>);
  m.def("ConnectedSubgraphs", &ConnectedSubgraphs<V1,E1,S1,D,P,P>);
  m.def("ConnectedSubgraphs", &ConnectedSubgraphs<V2,E2,S2,D,P,P>);
  
  using CycContainV1 = CondensedMolecularGraph::CycleVertContain;
  using CycContainE1 = CondensedMolecularGraph::CycleEdgeContain;
  using CycContainV2 = MolecularGraph::CycleVertContain;
  using CycContainE2 = MolecularGraph::CycleEdgeContain;
  m.def("CycleBasis", &CycleBasis<V1,E1,S1,D,P,P,CycContainV1>);
  m.def("CycleBasis", &CycleBasis<V1,E1,S1,D,P,P,CycContainE1>);
  m.def("CycleBasis", &CycleBasis<V2,E2,S2,D,P,P,CycContainV2>);
  m.def("CycleBasis", &CycleBasis<V2,E2,S2,D,P,P,CycContainE2>);
  m.def("AllCycles", &AllCycles<V1,E1,S1,D,P,P,CycContainE1>);
  m.def("AllCycles", &AllCycles<V2,E2,S2,D,P,P,CycContainE2>);
  
  m.def("ShortestPath", &ShortestPath<V1,E1,S1,D,P,P>);
  m.def("ShortestPath", &ShortestPath<V2,E2,S2,D,P,P>);
  m.def("AllSimplePaths", &AllSimplePaths<V1,E1,S1,D,P,P>);
  m.def("AllSimplePaths", &AllSimplePaths<V2,E2,S2,D,P,P>);
}
