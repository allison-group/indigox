#include <indigox/python/interface.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pyindigox, m) {
//  GenerateOpaqueContainers(m);
  // Base namespace
  GeneratePyMolecule(m);
  GeneratePyPeriodicTable(m);
  GeneratePyForcefield(m);
//  GeneratePyAthenaeum(m);
  // Graph namespace
  pybind11::module m_graph = m.def_submodule("graph");
  GeneratePyGraphs(m_graph);
//  GeneratePyElectronAssignmentGraph(m_graph);
//  GeneratePyCondensedMolecularGraph(m_graph);
  // Algorthm namespace
  pybind11::module m_algo = m.def_submodule("algorithm");
//  GeneratePyElectronAssigner(m_algo);
  
  // Minor utils things
//  py::enum_<indigox::utils::Option>(m, "Option")
//  .value("Yes", indigox::utils::Option::Yes)
//  .value("No", indigox::utils::Option::No)
//  .value("Auto", indigox::utils::Option::Auto)
//  .value("Default", indigox::utils::Option::Default)
//  .value("All", indigox::utils::Option::All)
//  .value("Some", indigox::utils::Option::Some)
//  .value("None", indigox::utils::Option::None)
//  ;
  
}


