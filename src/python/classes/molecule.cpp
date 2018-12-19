#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/molecular.hpp>

namespace py = pybind11;

void GeneratePyMolecule(pybind11::module& m) {
  using namespace indigox;
  py::return_value_policy Ref = py::return_value_policy::reference;
  // ===========================================================================
  // == Enum Types bindings ====================================================
  // ===========================================================================
  py::enum_<AtomStereo>(m, "AtomStereo")
  .value("Undefined", AtomStereo::UNDEFINED)
  .value("Achiral", AtomStereo::ACHIRAL)
  .value("R", AtomStereo::R)
  .value("S", AtomStereo::S)
  ;
  
  py::enum_<BondOrder>(m, "BondOrder")
  .value("Undefined", BondOrder::UNDEFINED)
  .value("Single", BondOrder::SINGLE)
  .value("Double", BondOrder::DOUBLE)
  .value("Triple", BondOrder::TRIPLE)
  .value("Quadruple", BondOrder::QUADRUPLE)
  .value("Aromatic", BondOrder::AROMATIC)
  .value("OneAndAHalf", BondOrder::ONEANDAHALF)
  .value("TwoAndAHalf", BondOrder::TWOANDAHALF)
  ;
  
  py::enum_<BondStereo>(m, "BondStereo")
  .value("Undefined", BondStereo::UNDEFINED)
  .value("None", BondStereo::NONE)
  .value("E", BondStereo::E)
  .value("Z", BondStereo::Z)
  ;
  
  // ===========================================================================
  // == Atom class bindings ====================================================
  // ===========================================================================
  py::class_<Atom>(m, "Atom")
  .def(py::init<>())
  .def("NumBonds", &Atom::NumBonds)
  .def("NumAngles", &Atom::NumAngles)
  .def("NumDihedrals", &Atom::NumDihedrals)
  .def("HasType", &Atom::HasType)
  .def("GetElement", &Atom::GetElement)
  .def("GetFormalCharge", &Atom::GetFormalCharge)
  .def("GetPartialCharge", &Atom::GetPartialCharge)
  .def("GetTag", &Atom::GetTag)
  .def("GetID", &Atom::GetID)
  .def("GetImplicitCount", &Atom::GetImplicitCount)
  .def("GetMolecule", &Atom::GetMolecule)
  .def("GetName", &Atom::GetName)
  .def("GetX", &Atom::GetX)
  .def("GetY", &Atom::GetY)
  .def("GetZ", &Atom::GetZ)
  .def("GetStereochemistry", &Atom::GetStereochemistry)
  .def("GetPosition", &Atom::GetPosition)
  .def("AddImplicitHydrogen", &Atom::AddImplicitHydrogen)
  .def("RemoveImplicitHydrogen", &Atom::RemoveImplicitHydrogen)
  .def("SetElement", py::overload_cast<const Element&>(&Atom::SetElement))
  .def("SetElement", py::overload_cast<std::string>(&Atom::SetElement))
  .def("SetElement", py::overload_cast<int32_t>(&Atom::SetElement))
  .def("SetFormalCharge", &Atom::SetFormalCharge)
  .def("SetPartialCharge", &Atom::SetPartialCharge)
  .def("SetImplicitCount", &Atom::SetImplicitCount)
  .def("SetTag", &Atom::SetTag)
  .def("SetName", &Atom::SetName)
  .def("SetPosition", &Atom::SetPosition)
  .def("SetStereochemistry", &Atom::SetStereochemistry)
  .def("GetBonds", &Atom::GetBonds, Ref)
  .def("GetAngles", &Atom::GetAngles, Ref)
  .def("GetDihedrals", &Atom::GetDihedrals, Ref)
  .def("GetIndex", &Atom::GetIndex)
  .def("GetType", &Atom::GetType)
  .def("SetType", &Atom::SetType)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &Atom::operator bool)
  .def("__str__", &outstream_operator<Atom>)
  .def("__repr__", &outstream_operator<Atom>)
  ;
  
  // ===========================================================================
  // == Bond class bindings ====================================================
  // ===========================================================================
  
  py::class_<Bond>(m, "Bond")
  .def(py::init<>())
  .def("GetTag", &Bond::GetTag)
  .def("GetID", &Bond::GetID)
  .def("GetMolecule", &Bond::GetMolecule)
  .def("GetOrder", &Bond::GetOrder)
  .def("GetStereochemistry", &Bond::GetStereochemistry)
  .def("SetTag", &Bond::SetTag)
  .def("SetOrder", &Bond::SetOrder)
  .def("GetAtoms", &Bond::GetAtoms)
  .def("NumAtoms", &Bond::NumAtoms)
  .def("GetIndex", &Bond::GetIndex)
  .def("GetType", &Bond::GetType)
  .def("SetType", &Bond::SetType)
  .def("HasType", &Bond::HasType)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &Bond::operator bool)
  .def("__str__", &outstream_operator<Bond>)
  .def("__repr__", &outstream_operator<Bond>)
  ;
  
  // ===========================================================================
  // == Angle class bindings ===================================================
  // ===========================================================================
  
  py::class_<Angle>(m, "Angle")
  .def(py::init<>())
  .def("HasType", &Angle::HasType)
  .def("NumAtoms", &Angle::NumAtoms)
  .def("GetTag", &Angle::GetTag)
  .def("GetID", &Angle::GetID)
  .def("GetMolecule", &Angle::GetMolecule)
  .def("GetAtoms", &Angle::GetAtoms)
  .def("GetIndex", &Angle::GetIndex)
  .def("GetType", &Angle::GetType)
  .def("SetType", &Angle::SetType)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &Angle::operator bool)
  .def("__str__", &outstream_operator<Angle>)
  .def("__repr__", &outstream_operator<Angle>)
  ;
  
  // ===========================================================================
  // == Dihedral class bindings ================================================
  // ===========================================================================
  
  py::class_<Dihedral>(m, "Dihedral")
  .def("GetTag", &Dihedral::GetTag)
  .def("GetID", &Dihedral::GetID)
  .def("GetMolecule", &Dihedral::GetMolecule)
  .def("GetAtoms", &Dihedral::GetAtoms)
  .def("NumAtoms", &Dihedral::NumAtoms)
  .def("SetTag", &Dihedral::SetTag)
  .def("GetIndex", &Dihedral::GetIndex)
  .def("NumTypes", &Dihedral::NumTypes)
  .def("GetTypes", &Dihedral::GetTypes, Ref)
  .def("HasType", &Dihedral::HasType)
  .def("SetTypes", &Dihedral::SetTypes)
  .def("AddType", &Dihedral::AddType)
  .def("RemoveType", &Dihedral::RemoveType)
  .def("GetPriority", &Dihedral::GetPriority)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &Dihedral::operator bool)
  .def("__str__", &outstream_operator<Dihedral>)
  .def("__repr__", &outstream_operator<Dihedral>)
  ;
  
  // ===========================================================================
  // == Molecule class bindings ================================================
  // ===========================================================================
  py::class_<Molecule>(m, "Molecule")
  .def(py::init<>())
  .def(py::init<std::string>())
  .def("HasAtom", &Molecule::HasAtom)
  .def("HasBond", py::overload_cast<const Bond&>(&Molecule::HasBond, py::const_))
  .def("HasBond", py::overload_cast<const Atom&, const Atom&>(&Molecule::HasBond, py::const_))
  .def("HasAngle", py::overload_cast<const Angle&>(&Molecule::HasAngle, py::const_))
  .def("HasAngle", py::overload_cast<const Atom&, const Atom&, const Atom&>(&Molecule::HasAngle))
  .def("HasDihedral", py::overload_cast<const Dihedral&>(&Molecule::HasDihedral, py::const_))
  .def("HasDihedral", py::overload_cast<const Atom&, const Atom&, const Atom&, const Atom&>(&Molecule::HasDihedral))
  .def("IsFrozen", &Molecule::IsFrozen)
  .def("NumAtoms", &Molecule::NumAtoms)
  .def("NumBonds", &Molecule::NumBonds)
  .def("NumAngles", &Molecule::NumAngles)
  .def("NumDihedrals", &Molecule::NumDihedrals)
  .def("GetAtom", &Molecule::GetAtom)
  .def("GetAtomID", &Molecule::GetAtomID)
  .def("GetAtomTag", &Molecule::GetAtomTag)
  .def("GetBond", py::overload_cast<uint32_t>(&Molecule::GetBond, py::const_))
  .def("GetBond", py::overload_cast<const Atom&, const Atom&>(&Molecule::GetBond, py::const_))
  .def("GetBondID", &Molecule::GetBondID)
  .def("GetBondTag", &Molecule::GetBondTag)
  .def("GetAngle", py::overload_cast<uint32_t>(&Molecule::GetAngle))
  .def("GetAngle", py::overload_cast<const Atom&, const Atom&, const Atom&>(&Molecule::GetAngle))
  .def("GetAngleID", &Molecule::GetAngleID)
  .def("GetAngleTag", &Molecule::GetAngleTag)
  .def("GetDihedral", py::overload_cast<uint32_t>(&Molecule::GetDihedral))
  .def("GetDihedral", py::overload_cast<const Atom&, const Atom&, const Atom&, const Atom&>(&Molecule::GetDihedral))
  .def("GetDihedralID", &Molecule::GetDihedralID)
  .def("GetDihedralTag", &Molecule::GetDihedralTag)
  .def("GetFormula", &Molecule::GetFormula)
  .def("GetGraph", &Molecule::GetGraph)
  .def("GetName", &Molecule::GetName)
  .def("GetMolecularCharge", &Molecule::GetMolecularCharge)
  .def("SetName", &Molecule::SetName)
  .def("SetMolecularCharge", &Molecule::SetMolecularCharge)
  
  .def("NewAtom", py::overload_cast<>(&Molecule::NewAtom))
  .def("NewAtom", py::overload_cast<const Element&>(&Molecule::NewAtom))
  .def("NewAtom", py::overload_cast<const Element&, double, double, double>(&Molecule::NewAtom))
  .def("NewBond", &Molecule::NewBond)
  .def("NewDihedral", py::overload_cast<const Atom&, const Atom&, const Atom&, const Atom&>(&Molecule::NewDihedral))
  .def("RemoveAtom", &Molecule::RemoveAtom)
  .def("RemoveBond", py::overload_cast<const Bond&>(&Molecule::RemoveBond))
  .def("RemoveBond", py::overload_cast<const Atom&, const Atom&>(&Molecule::RemoveBond))
  .def("PerceiveAngles", &Molecule::PerceiveAngles)
  .def("PerceiveDihedrals", &Molecule::PerceiveDihedrals)
  .def("ReserveAtoms", &Molecule::ReserveAtoms)
  .def("ReserveBonds", &Molecule::ReserveBonds)
  .def("GetAtoms", &Molecule::GetAtoms, Ref)
  .def("GetBonds", &Molecule::GetBonds, Ref)
  .def("GetAngles", &Molecule::GetAngles, Ref)
  .def("GetDihedrals", &Molecule::GetDihedrals, Ref)
  .def("GetForcefield", &Molecule::GetForcefield)
  .def("SetForcefield", &Molecule::SetForcefield)
  .def("ResetForcefield", &Molecule::ResetForcefield)
  .def("HasForcefield", &Molecule::HasForcefield)
  .def("GetCurrentState", &Molecule::GetCurrentState)
  .def("ModificationMade", &Molecule::ModificationMade)
  .def("FreezeModifications", &Molecule::FreezeModifications)
  ;
  
  // ===========================================================================
  // == Module function bindings ===============================================
  // ===========================================================================
  m.def("SaveMolecule", &SaveMolecule);
  m.def("LoadMolecule", &LoadMolecule);
  m.def("Benzene", []() -> Molecule {
    Molecule mol("Benzene");
    Atom C1 = mol.NewAtom(GetPeriodicTable().GetElement("C"));
    Atom C2 = mol.NewAtom(GetPeriodicTable().GetElement("C"));
    Atom C3 = mol.NewAtom(GetPeriodicTable().GetElement("C"));
    Atom C4 = mol.NewAtom(GetPeriodicTable().GetElement("C"));
    Atom C5 = mol.NewAtom(GetPeriodicTable().GetElement("C"));
    Atom C6 = mol.NewAtom(GetPeriodicTable().GetElement("C"));
    Atom H1 = mol.NewAtom(GetPeriodicTable().GetElement("H"));
    Atom Cl = mol.NewAtom(GetPeriodicTable().GetElement("Cl"));
    Atom F1 = mol.NewAtom(GetPeriodicTable().GetElement("F"));
    Atom Br = mol.NewAtom(GetPeriodicTable().GetElement("Br"));
    Atom I1 = mol.NewAtom(GetPeriodicTable().GetElement("I"));
    Atom H2 = mol.NewAtom(GetPeriodicTable().GetElement("H"));
    mol.NewBond(C1, C2).SetOrder(BondOrder::AROMATIC);
    mol.NewBond(C2, C3).SetOrder(BondOrder::AROMATIC);
    mol.NewBond(C3, C4).SetOrder(BondOrder::AROMATIC);
    mol.NewBond(C4, C5).SetOrder(BondOrder::AROMATIC);
    mol.NewBond(C5, C6).SetOrder(BondOrder::AROMATIC);
    mol.NewBond(C6, C1).SetOrder(BondOrder::AROMATIC);
    mol.NewBond(C1, H1).SetOrder(BondOrder::SINGLE);
    mol.NewBond(C2, Cl).SetOrder(BondOrder::SINGLE);
    mol.NewBond(C3, F1).SetOrder(BondOrder::SINGLE);
    mol.NewBond(C4, Br).SetOrder(BondOrder::SINGLE);
    mol.NewBond(C5, I1).SetOrder(BondOrder::SINGLE);
    mol.NewBond(C6, H2).SetOrder(BondOrder::SINGLE);
    mol.FreezeModifications();
    return mol;
  });
  
  // Container bindings
  py::bind_vector<std::vector<Atom>>(m, "VecAtom");
  py::bind_vector<std::vector<Bond>>(m, "VecBond");
  py::bind_vector<std::vector<Angle>>(m, "VecAngle");
  py::bind_vector<std::vector<Dihedral>>(m, "VecDihedral");

}
