#include <indigox/python/interface.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pyindigox, m) {
  GenerateOptions(m);
  // Base namespace
  GeneratePyAngle(m);
  GeneratePyAtom(m);
  GeneratePyBond(m);
  GeneratePyDihedral(m);
  GeneratePyElement(m);
  GeneratePyMolecule(m);
  GeneratePyPeriodicTable(m);
  // Graph namespace
  pybind11::module m_graph = m.def_submodule("graph");
  GeneratePyMolecularGraph(m_graph);
  GeneratePyElectronAssignmentGraph(m_graph);
}


