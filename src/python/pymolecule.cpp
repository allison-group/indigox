#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include <indigox/python/interface.hpp>

#include <indigox/classes/molecule.hpp>

namespace py = pybind11;

void GeneratePyMolecule(py::module& m) {
  using namespace indigox;
  
  m.def("CreateMolecule", &CreateMolecule);
  
  auto angle_get = [](const Molecule& mol) {
    auto iters = mol->GetAngleIters();
    return py::make_iterator(iters.first, iters.second);
  };
  
  auto atom_get = [](const Molecule& mol) {
    auto iters = mol->GetAtomIters();
    return py::make_iterator(iters.first, iters.second);
  };
  
  auto bond_get = [](const Molecule& mol) {
    auto iters = mol->GetBondIters();
    return py::make_iterator(iters.first, iters.second);
  };
  
  auto dihedral_get = [](const Molecule& mol) {
    auto iters = mol->GetDihedralIters();
    return py::make_iterator(iters.first, iters.second);
  };
  
  using crAtom = const Atom&;
  
  py::class_<IXMolecule, Molecule>(m, "Molecule")
  // Hidden methods
  //  No constructor
  .def("__repr__", [](Molecule mol) {
    std::stringstream ss; ss << mol; return ss.str();
  })
  // Getters
  .def("GetAngle", py::overload_cast<size_t>(&IXMolecule::GetAngle, py::const_))
  .def("GetAngle", py::overload_cast<crAtom, crAtom, crAtom>(&IXMolecule::GetAngle))
  .def("GetAngleID", &IXMolecule::GetAngleID)
  .def("GetAngles", angle_get, py::keep_alive<0, 1>())
  .def("GetAngleTag", &IXMolecule::GetAngleTag)
  .def("GetAtom", &IXMolecule::GetAtom)
  .def("GetAtomID", &IXMolecule::GetAtomID)
  .def("GetAtoms", atom_get, py::keep_alive<0, 1>())
  .def("GetAtomTag", &IXMolecule::GetAtomTag)
  .def("GetBond", py::overload_cast<size_t>(&IXMolecule::GetBond, py::const_))
  .def("GetBond", py::overload_cast<crAtom, crAtom>(&IXMolecule::GetBond, py::const_))
  .def("GetBondID", &IXMolecule::GetBondID)
  .def("GetBonds", bond_get, py::keep_alive<0, 1>())
  .def("GetBondTag", &IXMolecule::GetBondTag)
  .def("GetDihedral", py::overload_cast<size_t>(&IXMolecule::GetDihedral, py::const_))
  .def("GetDihedral", py::overload_cast<crAtom, crAtom, crAtom, crAtom>(&IXMolecule::GetDihedral))
  .def("GetDihedralID", &IXMolecule::GetDihedralID)
  .def("GetDihedrals", dihedral_get, py::keep_alive<0, 1>())
  .def("GetDihedralTag", &IXMolecule::GetDihedralTag)
  .def("GetFormula", &IXMolecule::GetFormula)
  .def("GetGraph", &IXMolecule::GetGraph)
  .def("GetMolecularCharge", &IXMolecule::GetMolecularCharge)
  .def("GetName", &IXMolecule::GetName)
  .def("GetUniqueID", &IXMolecule::GetUniqueID)
  .def("NumAngles", &IXMolecule::NumAngles)
  .def("NumAtoms", &IXMolecule::NumAtoms)
  .def("NumBonds", &IXMolecule::NumBonds)
  .def("NumDihedrals", &IXMolecule::NumDihedrals)
  // Checkers
  .def("HasAngle", py::overload_cast<const Angle&>(&IXMolecule::HasAngle, py::const_))
  .def("HasAngle", py::overload_cast<crAtom, crAtom, crAtom>(&IXMolecule::HasAngle))
  .def("HasAtom", &IXMolecule::HasAtom)
  .def("HasBond", py::overload_cast<const Bond&>(&IXMolecule::HasBond, py::const_))
  .def("HasBond", py::overload_cast<crAtom, crAtom>(&IXMolecule::HasBond, py::const_))
  .def("HasDihedral", py::overload_cast<const Dihedral&>(&IXMolecule::HasDihedral, py::const_))
  .def("HasDihedral", py::overload_cast<crAtom, crAtom, crAtom, crAtom>(&IXMolecule::HasDihedral))
  // Data modification
  .def("SetMolecularCharge", &IXMolecule::SetMolecularCharge)
  .def("SetName", &IXMolecule::SetName)
  .def("SetPropertyModified", &IXMolecule::SetPropertyModified)
  // Molecule modification
  .def("NewAtom", py::overload_cast<>(&IXMolecule::NewAtom))
  .def("NewAtom", py::overload_cast<Element>(&IXMolecule::NewAtom))
  .def("NewAtom", py::overload_cast<std::string>(&IXMolecule::NewAtom))
  .def("NewAtom", py::overload_cast<std::string, Element>(&IXMolecule::NewAtom))
  .def("NewBond", &IXMolecule::NewBond)
  .def("PerceiveAngles", &IXMolecule::PerceiveAngles)
  .def("PerceiveDihedrals", &IXMolecule::PerceiveDihedrals)
  .def("RemoveAtom", &IXMolecule::RemoveAtom)
  .def("RemoveBond", py::overload_cast<Bond>(&IXMolecule::RemoveBond))
  .def("RemoveBond", py::overload_cast<Atom, Atom>(&IXMolecule::RemoveBond))
  // Other
  .def("ReserveAtoms", &IXMolecule::ReserveAtoms)
  .def("ReserveBonds", &IXMolecule::ReserveBonds)
  ;
  
  py::enum_<MolProperty>(m, "MolProperty")
  .value("ATOM_ELEMENTS", MolProperty::ATOM_ELEMENTS)
  .value("CONNECTIVITY", MolProperty::CONNECTIVITY)
  .value("ELECTRON_COUNT", MolProperty::ELECTRON_COUNT)
  ;
}
