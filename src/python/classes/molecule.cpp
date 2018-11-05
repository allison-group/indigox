#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>

namespace py = pybind11;

void GeneratePyMolecule(pybind11::module& m) {
  using namespace indigox;
  py::return_value_policy Ref = py::return_value_policy::reference;
  // ===========================================================================
  // == Enum Types bindings ====================================================
  // ===========================================================================
  py::enum_<AtomStereo>(m, "AtomStereo")
  .value("UNDEFINED", AtomStereo::UNDEFINED)
  .value("ACHIRAL", AtomStereo::ACHIRAL)
  .value("R", AtomStereo::R)
  .value("S", AtomStereo::S)
  ;
  
  py::enum_<BondOrder>(m, "BondOrder")
  .value("UNDEFINED", BondOrder::UNDEFINED)
  .value("SINGLE", BondOrder::SINGLE)
  .value("DOUBLE", BondOrder::DOUBLE)
  .value("TRIPLE", BondOrder::TRIPLE)
  .value("QUADRUPLE", BondOrder::QUADRUPLE)
  .value("AROMATIC", BondOrder::AROMATIC)
  .value("ONEANDAHALF", BondOrder::ONEANDAHALF)
  .value("TWOANDAHALF", BondOrder::TWOANDAHALF)
  ;
  
  py::enum_<BondStereo>(m, "BondStereo")
  .value("UNDEFINED", BondStereo::UNDEFINED)
  .value("NONE", BondStereo::NONE)
  .value("E", BondStereo::E)
  .value("Z", BondStereo::Z)
  ;
  
  // ===========================================================================
  // == Atom class bindings ====================================================
  // ===========================================================================
  py::class_<Atom, sAtom>(m, "Atom")
  .def("GetElement", &Atom::GetElement)
  .def("GetFormalCharge", &Atom::GetFormalCharge)
  .def("GetPartialCharge", &Atom::GetPartialCharge)
  .def("GetTag", &Atom::GetTag)
  .def("GetImplicitCount", &Atom::GetImplicitCount)
  .def("AddImplicitHydrogen", &Atom::AddImplicitHydrogen)
  .def("RemoveImplicitHydrogen", &Atom::RemoveImplicitHydrogen)
  .def("GetMolecule", &Atom::GetMolecule, Ref)
  .def("GetName", &Atom::GetName)
  .def("GetX", &Atom::GetX)
  .def("GetY", &Atom::GetY)
  .def("GetZ", &Atom::GetZ)
  .def("GetVector", &Atom::GetVector)
  .def("ToString", &Atom::ToString)
  .def("SetElement", py::overload_cast<const Element&>(&Atom::SetElement))
  .def("SetElement", py::overload_cast<std::string>(&Atom::SetElement))
  .def("SetElement", py::overload_cast<unsigned int>(&Atom::SetElement))
  .def("SetFormalCharge", &Atom::SetFormalCharge)
  .def("SetPartialCharge", &Atom::SetPartialCharge)
  .def("SetImplicitCount", &Atom::SetImplicitCount)
  .def("SetTag", &Atom::SetTag)
  .def("SetName", &Atom::SetName)
  .def("SetX", &Atom::SetX)
  .def("SetY", &Atom::SetY)
  .def("SetZ", &Atom::SetZ)
  .def("SetPosition", &Atom::SetPosition)
  .def("SetStereochemistry", &Atom::SetStereochemistry)
  .def("SetAromaticity", &Atom::SetAromaticity)
  .def("GetStereochemistry", &Atom::GetStereochemistry)
  .def("GetAromaticity", &Atom::GetAromaticity)
  .def("GetBonds", &Atom::GetBonds, Ref)
  .def("NumBonds", &Atom::NumBonds)
  .def("NumAngles", &Atom::NumAngles)
  .def("NumDihedrals", &Atom::NumDihedrals)
  .def("GetIndex", &Atom::GetIndex)
  .def("GetType", &Atom::GetType)
  .def("HasType", &Atom::HasType)
  .def("SetType", &Atom::SetType)
  .def("__str__", &outstream_operator<Atom>)
  .def("__repr__", &outstream_operator<Atom>)
  
  .def("GetUniqueID", &Atom::GetUniqueID)
  ;
  
  // ===========================================================================
  // == Bond class bindings ====================================================
  // ===========================================================================
  py::class_<Bond, sBond>(m, "Bond")
  .def("GetTag", &Bond::GetTag)
  .def("GetMolecule", &Bond::GetMolecule, Ref)
  .def("GetOrder", &Bond::GetOrder)
  .def("GetSourceAtom", &Bond::GetSourceAtom, Ref)
  .def("GetTargetAtom", &Bond::GetTargetAtom, Ref)
  .def("GetAromaticity", &Bond::GetAromaticity)
  .def("GetStereochemistry", &Bond::GetStereochemistry)
  .def("ToString", &Bond::ToString)
  .def("SetTag", &Bond::SetTag)
  .def("SetOrder", &Bond::SetOrder)
  .def("SwapOrder", &Bond::SwapOrder)
  .def("SetAromaticity", &Bond::SetAromaticity)
  .def("SetStereochemistry", &Bond::SetStereochemistry)
  .def("GetAtoms", &Bond::GetAtoms, Ref)
  .def("NumAtoms", &Bond::NumAtoms)
  .def("GetIndex", &Bond::GetIndex)
  .def("GetType", &Bond::GetType)
  .def("SetType", &Bond::SetType)
  .def("HasType", &Bond::HasType)
  .def("__str__", &outstream_operator<Bond>)
  .def("__repr__", &outstream_operator<Bond>)
  
  .def("GetUniqueID", &Bond::GetUniqueID)
  ;
  
  // ===========================================================================
  // == Angle class bindings ===================================================
  // ===========================================================================
  auto get_atoms_ang = [](const Angle& ang) {
    stdx::triple<Atom&, Atom&, Atom&> atms = ang.GetAtoms();
    return std::make_tuple(atms.first, atms.second, atms.third);
  };
  
  py::class_<Angle, sAngle>(m, "Angle")
  .def("GetTag", &Angle::GetTag)
  .def("GetMolecule", &Angle::GetMolecule, Ref)
  .def("ToString", &Angle::ToString)
  .def("SetTag", &Angle::SetTag)
  .def("SwapOrder", &Angle::SwapOrder)
  .def("GetAtoms", get_atoms_ang, Ref)
  .def("NumAtoms", &Angle::NumAtoms)
  .def("GetIndex", &Angle::GetIndex)
  .def("GetType", &Angle::GetType)
  .def("SetType", &Angle::SetType)
  .def("HasType", &Angle::HasType)
  .def("__str__", &outstream_operator<Angle>)
  .def("__repr__", &outstream_operator<Angle>)
  
  .def("GetUniqueID", &Angle::GetUniqueID)
  ;
  
  // ===========================================================================
  // == Dihedral class bindings ================================================
  // ===========================================================================
  auto get_atoms_dhd = [](const Dihedral& dhd) {
    stdx::quad<Atom&, Atom&, Atom&, Atom&> atms = dhd.GetAtoms();
    return std::make_tuple(atms.first, atms.second, atms.third, atms.fourth);
  };
  
  py::class_<Dihedral, sDihedral>(m, "Dihedral")
  .def("GetTag", &Dihedral::GetTag)
  .def("GetMolecule", &Dihedral::GetMolecule, Ref)
  .def("GetAtoms", get_atoms_dhd, Ref)
  .def("NumAtoms", &Dihedral::NumAtoms)
  .def("SwapOrder", &Dihedral::SwapOrder)
  .def("ToString", &Dihedral::ToString)
  .def("SetTag", &Dihedral::SetTag)
  .def("GetIndex", &Dihedral::GetIndex)
  .def("GetType", &Dihedral::GetType)
  .def("NumTypes", &Dihedral::NumTypes)
  .def("GetTypes", &Dihedral::GetTypes, Ref)
  .def("HasType", &Dihedral::HasType)
  .def("AddType", &Dihedral::AddType)
  .def("RemoveType", &Dihedral::RemoveType)
  .def("__str__", &outstream_operator<Dihedral>)
  .def("__repr__", &outstream_operator<Dihedral>)
  
  .def("GetUniqueID", &Dihedral::GetUniqueID)
  ;
  
  // ===========================================================================
  // == Molecule class bindings ================================================
  // ===========================================================================
  py::class_<Molecule, sMolecule>(m, "Molecule")
  .def("__init__", &CreateMolecule)
  .def("GetAngle", py::overload_cast<size_t>(&Molecule::GetAngle, py::const_), Ref)
  .def("GetDihedral", py::overload_cast<size_t>(&Molecule::GetDihedral, py::const_), Ref)
  .def("GetAngle", py::overload_cast<Atom&, Atom&, Atom&>(&Molecule::GetAngle), Ref)
  .def("GetDihedral", py::overload_cast<Atom&, Atom&, Atom&, Atom&>(&Molecule::GetDihedral), Ref)
  .def("GetAngleTag", &Molecule::GetAngleTag, Ref)
  .def("GetDihedralTag", &Molecule::GetDihedralTag, Ref)
  .def("GetAngleID", &Molecule::GetAngleID, Ref)
  .def("GetDihedralID", &Molecule::GetDihedralID, Ref)
  .def("GetAtom", &Molecule::GetAtom, Ref)
  .def("GetAtomTag", &Molecule::GetAtomTag, Ref)
  .def("GetAtomID", &Molecule::GetAtomID, Ref)
  .def("GetBond", py::overload_cast<size_t>(&Molecule::GetBond, py::const_), Ref)
  .def("GetBond", py::overload_cast<Atom&, Atom&>(&Molecule::GetBond, py::const_), Ref)
  .def("GetBondTag", &Molecule::GetBondTag, Ref)
  .def("GetBondID", &Molecule::GetBondID, Ref)
  .def("GetFormula", &Molecule::GetFormula)
  .def("GetGraph", &Molecule::GetGraph, Ref)
  .def("GetName", &Molecule::GetName)
  .def("GetMolecularCharge", &Molecule::GetMolecularCharge)
  .def("NumAtoms", &Molecule::NumAtoms)
  .def("NumBonds", &Molecule::NumBonds)
  .def("NumAngles", &Molecule::NumAngles)
  .def("NumDihedrals", &Molecule::NumDihedrals)
  .def("SetName", &Molecule::SetName)
  .def("SetMolecularCharge", &Molecule::SetMolecularCharge)
  .def("HasAtom", &Molecule::HasAtom)
  .def("HasBond", py::overload_cast<Bond&>(&Molecule::HasBond, py::const_))
  .def("HasBond", py::overload_cast<Atom&, Atom&>(&Molecule::HasBond, py::const_))
  .def("HasAngle", py::overload_cast<Angle&>(&Molecule::HasAngle, py::const_))
  .def("HasAngle", py::overload_cast<Atom&, Atom&, Atom&>(&Molecule::HasAngle))
  .def("HasDihedral", py::overload_cast<Dihedral&>(&Molecule::HasDihedral, py::const_))
  .def("HasDihedral", py::overload_cast<Atom&, Atom&, Atom&, Atom&>(&Molecule::HasDihedral))
  .def("NewAtom", py::overload_cast<>(&Molecule::NewAtom), Ref)
  .def("NewAtom", py::overload_cast<const Element&>(&Molecule::NewAtom), Ref)
  .def("NewAtom", py::overload_cast<std::string>(&Molecule::NewAtom), Ref)
  .def("NewAtom", py::overload_cast<std::string, const Element&>(&Molecule::NewAtom), Ref)
  .def("NewBond", &Molecule::NewBond, Ref)
  .def("RemoveAtom", &Molecule::RemoveAtom)
  .def("RemoveBond", py::overload_cast<Bond&>(&Molecule::RemoveBond))
  .def("RemoveBond", py::overload_cast<Atom&, Atom&>(&Molecule::RemoveBond))
  .def("NewDihedral", &Molecule::NewDihedral, Ref)
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
  .def("HasForcefield", &Molecule::HasForcefield)
  
  .def("GetUniqueID", &Molecule::GetUniqueID)
  
  .def("GetCurrentState", &Molecule::GetCurrentState)
  .def("ModificationMade", &Molecule::ModificationMade)
  .def("FreezeModifications", &Molecule::FreezeModifications)
  .def("IsFrozen", &Molecule::IsFrozen)
  ;
  
  // ===========================================================================
  // == Module function bindings ===============================================
  // ===========================================================================
  m.def("CreateMolecule", &CreateMolecule);
  m.def("Benzene", []() -> sMolecule {
    sMolecule mol = CreateMolecule();
    Atom& C1 = mol->NewAtom("C1", GetPeriodicTable().GetElement("C"));
    Atom& C2 = mol->NewAtom("C2", GetPeriodicTable().GetElement("C"));
    Atom& C3 = mol->NewAtom("C3", GetPeriodicTable().GetElement("C"));
    Atom& C4 = mol->NewAtom("C4", GetPeriodicTable().GetElement("C"));
    Atom& C5 = mol->NewAtom("C5", GetPeriodicTable().GetElement("C"));
    Atom& C6 = mol->NewAtom("C6", GetPeriodicTable().GetElement("C"));
    Atom& H1 = mol->NewAtom("H1", GetPeriodicTable().GetElement("H"));
    Atom& Cl = mol->NewAtom("Cl", GetPeriodicTable().GetElement("Cl"));
    Atom& F1 = mol->NewAtom("F1", GetPeriodicTable().GetElement("F"));
    Atom& Br = mol->NewAtom("Br", GetPeriodicTable().GetElement("Br"));
    Atom& I1 = mol->NewAtom("I1", GetPeriodicTable().GetElement("I"));
    Atom& H2 = mol->NewAtom("H2", GetPeriodicTable().GetElement("H"));
    mol->NewBond(C1, C2).SetOrder(BondOrder::AROMATIC);
    mol->NewBond(C2, C3).SetOrder(BondOrder::AROMATIC);
    mol->NewBond(C3, C4).SetOrder(BondOrder::AROMATIC);
    mol->NewBond(C4, C5).SetOrder(BondOrder::AROMATIC);
    mol->NewBond(C5, C6).SetOrder(BondOrder::AROMATIC);
    mol->NewBond(C6, C1).SetOrder(BondOrder::AROMATIC);
    mol->NewBond(C1, H1).SetOrder(BondOrder::SINGLE);
    mol->NewBond(C2, Cl).SetOrder(BondOrder::SINGLE);
    mol->NewBond(C3, F1).SetOrder(BondOrder::SINGLE);
    mol->NewBond(C4, Br).SetOrder(BondOrder::SINGLE);
    mol->NewBond(C5, I1).SetOrder(BondOrder::SINGLE);
    mol->NewBond(C6, H2).SetOrder(BondOrder::SINGLE);
    mol->FreezeModifications();
    return mol;
  });
}
