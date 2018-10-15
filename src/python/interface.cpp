#include <indigox/python/interface.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/periodictable.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pyindigox, m) {
  GenerateOpaqueContainers(m);
  // Base namespace
  GeneratePyAngle(m);
  GeneratePyAtom(m);
  GeneratePyBond(m);
  GeneratePyDihedral(m);
  GeneratePyElement(m);
  GeneratePyMolecule(m);
  GeneratePyPeriodicTable(m);
  GeneratePyForcefield(m);
  GeneratePyAthenaeum(m);
  // Graph namespace
  pybind11::module m_graph = m.def_submodule("graph");
  GeneratePyMolecularGraph(m_graph);
  GeneratePyElectronAssignmentGraph(m_graph);
  GeneratePyCondensedMolecularGraph(m_graph);
  // Algorthm namespace
  pybind11::module m_algo = m.def_submodule("algorithm");
  GeneratePyElectronAssigner(m_algo);
  
  // Minor utils things
  py::enum_<indigox::utils::Option>(m, "Option")
  .value("Yes", indigox::utils::Option::Yes)
  .value("No", indigox::utils::Option::No)
  .value("Auto", indigox::utils::Option::Auto)
  .value("Default", indigox::utils::Option::Default)
  .value("All", indigox::utils::Option::All)
  .value("Some", indigox::utils::Option::Some)
  .value("None", indigox::utils::Option::None)
  ;
  
  m.def("Benzene", []() -> indigox::Molecule {
    using namespace indigox;
    Molecule mol = CreateMolecule();
    Atom C1 = mol->NewAtom("C1", GetPeriodicTable().GetElement("C"));
    Atom C2 = mol->NewAtom("C2", GetPeriodicTable().GetElement("C"));
    Atom C3 = mol->NewAtom("C3", GetPeriodicTable().GetElement("C"));
    Atom C4 = mol->NewAtom("C4", GetPeriodicTable().GetElement("C"));
    Atom C5 = mol->NewAtom("C5", GetPeriodicTable().GetElement("C"));
    Atom C6 = mol->NewAtom("C6", GetPeriodicTable().GetElement("C"));
    Atom H1 = mol->NewAtom("H1", GetPeriodicTable().GetElement("H"));
    Atom Cl = mol->NewAtom("Cl", GetPeriodicTable().GetElement("Cl"));
    Atom F1 = mol->NewAtom("F1", GetPeriodicTable().GetElement("F"));
    Atom Br = mol->NewAtom("Br", GetPeriodicTable().GetElement("Br"));
    Atom I1 = mol->NewAtom("I1", GetPeriodicTable().GetElement("I"));
    Atom H2 = mol->NewAtom("H2", GetPeriodicTable().GetElement("H"));
    mol->NewBond(C1, C2)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C2, C3)->SetOrder(BondOrder::DOUBLE);
    mol->NewBond(C3, C4)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C4, C5)->SetOrder(BondOrder::DOUBLE);
    mol->NewBond(C5, C6)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C6, C1)->SetOrder(BondOrder::DOUBLE);
    mol->NewBond(C1, H1)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C2, Cl)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C3, F1)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C4, Br)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C5, I1)->SetOrder(BondOrder::SINGLE);
    mol->NewBond(C6, H2)->SetOrder(BondOrder::SINGLE);
    return mol;
  });
}


