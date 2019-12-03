#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/residue.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/python/interface.hpp>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace special {
  indigox::Molecule Methionine() {
    using namespace indigox;
    
    PeriodicTable pt = GetPeriodicTable();
    Molecule mol("Methionine");
    
    Atom N = mol.NewAtom(pt["N"]); N.SetFormalCharge(1); N.SetName("N");
    Atom H1 = mol.NewAtom(pt["H"]); H1.SetFormalCharge(0); H1.SetName("H1");
    Atom H2 = mol.NewAtom(pt["H"]); H2.SetFormalCharge(0); H2.SetName("H2");
    Atom H3 = mol.NewAtom(pt["H"]); H3.SetFormalCharge(0); H3.SetName("H3");
    Atom CA = mol.NewAtom(pt["C"]); CA.SetFormalCharge(0); CA.SetName("CA");
    Atom HA1 = mol.NewAtom(pt["H"]); HA1.SetFormalCharge(0); HA1.SetName("HA1");
    Atom C = mol.NewAtom(pt["C"]); C.SetFormalCharge(0); C.SetName("C");
    Atom O1 = mol.NewAtom(pt["O"]); O1.SetFormalCharge(0); O1.SetName("O1");
    Atom O2 = mol.NewAtom(pt["O"]); O2.SetFormalCharge(-1); O2.SetName("O2");
    Atom CB = mol.NewAtom(pt["C"]); CB.SetFormalCharge(0); CB.SetName("CB");
    Atom HB1 = mol.NewAtom(pt["H"]); HB1.SetFormalCharge(0); HB1.SetName("HB2");
    Atom HB2 = mol.NewAtom(pt["H"]); HB2.SetFormalCharge(0); HB2.SetName("HB1");
    Atom CG = mol.NewAtom(pt["C"]); CG.SetFormalCharge(0); CG.SetName("CG");
    Atom HG1 = mol.NewAtom(pt["H"]); HG1.SetFormalCharge(0); HG1.SetName("HG2");
    Atom HG2 = mol.NewAtom(pt["H"]); HG2.SetFormalCharge(0); HG2.SetName("HG1");
    Atom SD = mol.NewAtom(pt["S"]); SD.SetFormalCharge(0); SD.SetName("SD");
    Atom CE = mol.NewAtom(pt["C"]); CE.SetFormalCharge(0); CE.SetName("CE");
    Atom HE1 = mol.NewAtom(pt["H"]); HE1.SetFormalCharge(0); HE1.SetName("HE2");
    Atom HE2 = mol.NewAtom(pt["H"]); HE2.SetFormalCharge(0); HE2.SetName("HE1");
    Atom HE3 = mol.NewAtom(pt["H"]); HE3.SetFormalCharge(0); HE3.SetName("HE3");
    
    mol.NewBond(N, H1);
    mol.NewBond(N, H2);
    mol.NewBond(N, H3);
    mol.NewBond(N, CA);
    mol.NewBond(CA, HA1);
    mol.NewBond(CA, C);
    mol.NewBond(CA, CB);
    mol.NewBond(C, O1).SetOrder(BondOrder::DOUBLE);
    mol.NewBond(C, O2);
    mol.NewBond(CB, HB1);
    mol.NewBond(CB, HB2);
    mol.NewBond(CB, CG);
    mol.NewBond(CG, HG1);
    mol.NewBond(CG, HG2);
    mol.NewBond(CG, SD);
    mol.NewBond(SD, CE);
    mol.NewBond(CE, HE1);
    mol.NewBond(CE, HE2);
    mol.NewBond(CE, HE3);
    
    return mol;
  }
  
  indigox::Molecule ProtonatedLysine() {
    using namespace indigox;
    PeriodicTable pt = GetPeriodicTable();
    Molecule mol("Lysine");
    
    Atom N = mol.NewAtom(pt["N"]); N.SetFormalCharge(1); N.SetName("N");
    Atom H1 = mol.NewAtom(pt["H"]); H1.SetFormalCharge(0); H1.SetName("H1");
    Atom H2 = mol.NewAtom(pt["H"]); H2.SetFormalCharge(0); H2.SetName("H2");
    Atom H3 = mol.NewAtom(pt["H"]); H3.SetFormalCharge(0); H3.SetName("H3");
    Atom CA = mol.NewAtom(pt["C"]); CA.SetFormalCharge(0); CA.SetName("CA");
    Atom HA1 = mol.NewAtom(pt["H"]); HA1.SetFormalCharge(0); HA1.SetName("HA1");
    Atom C = mol.NewAtom(pt["C"]); C.SetFormalCharge(0); C.SetName("C");
    Atom O1 = mol.NewAtom(pt["O"]); O1.SetFormalCharge(0); O1.SetName("O1");
    Atom O2 = mol.NewAtom(pt["O"]); O2.SetFormalCharge(-1); O2.SetName("O2");
    Atom CB = mol.NewAtom(pt["C"]); CB.SetFormalCharge(0); CB.SetName("CB");
    Atom HB1 = mol.NewAtom(pt["H"]); HB1.SetFormalCharge(0); HB1.SetName("HB2");
    Atom HB2 = mol.NewAtom(pt["H"]); HB2.SetFormalCharge(0); HB2.SetName("HB1");
    Atom CG = mol.NewAtom(pt["C"]); CG.SetFormalCharge(0); CG.SetName("CG");
    Atom HG1 = mol.NewAtom(pt["H"]); HG1.SetFormalCharge(0); HG1.SetName("HG2");
    Atom HG2 = mol.NewAtom(pt["H"]); HG2.SetFormalCharge(0); HG2.SetName("HG1");
    Atom CD = mol.NewAtom(pt["C"]); CD.SetFormalCharge(0); CD.SetName("CD");
    Atom HD1 = mol.NewAtom(pt["H"]); HD1.SetFormalCharge(0); HD1.SetName("HD1");
    Atom HD2 = mol.NewAtom(pt["H"]); HD2.SetFormalCharge(0); HD2.SetName("HD2");
    Atom CE = mol.NewAtom(pt["C"]); CE.SetFormalCharge(0); CE.SetName("CE");
    Atom HE1 = mol.NewAtom(pt["H"]); HE1.SetFormalCharge(0); HE1.SetName("HE2");
    Atom HE2 = mol.NewAtom(pt["H"]); HE2.SetFormalCharge(0); HE2.SetName("HE1");
    Atom NZ = mol.NewAtom(pt["N"]); NZ.SetFormalCharge(1); NZ.SetName("NZ");
    Atom HZ1 = mol.NewAtom(pt["H"]); HZ1.SetFormalCharge(0); HZ1.SetName("HZ1");
    Atom HZ2 = mol.NewAtom(pt["H"]); HZ2.SetFormalCharge(0); HZ2.SetName("HZ2");
    Atom HZ3 = mol.NewAtom(pt["H"]); HZ3.SetFormalCharge(0); HZ3.SetName("HZ3");
    
    mol.NewBond(N, H1);
    mol.NewBond(N, H2);
    mol.NewBond(N, H3);
    mol.NewBond(N, CA);
    mol.NewBond(CA, HA1);
    mol.NewBond(CA, C);
    mol.NewBond(CA, CB);
    mol.NewBond(C, O1).SetOrder(BondOrder::DOUBLE);
    mol.NewBond(C, O2);
    mol.NewBond(CB, HB1);
    mol.NewBond(CB, HB2);
    mol.NewBond(CB, CG);
    mol.NewBond(CG, HG1);
    mol.NewBond(CG, HG2);
    mol.NewBond(CG, CD);
    mol.NewBond(CD, HD1);
    mol.NewBond(CD, HD2);
    mol.NewBond(CD, CE);
    mol.NewBond(CE, HE1);
    mol.NewBond(CE, HE2);
    mol.NewBond(CE, NZ);
    mol.NewBond(NZ, HZ1);
    mol.NewBond(NZ, HZ2);
    mol.NewBond(NZ, HZ3);
    
    return mol;
  }
}

void GeneratePyMolecule(pybind11::module &m) {
  using namespace indigox;
  py::return_value_policy Ref = py::return_value_policy::reference;
  // ===========================================================================
  // == Enum Types bindings ====================================================
  // ===========================================================================
  py::enum_<AtomStereo>(m, "AtomStereo")
      .value("Undefined", AtomStereo::UNDEFINED)
      .value("Achiral", AtomStereo::ACHIRAL)
      .value("R", AtomStereo::R)
      .value("S", AtomStereo::S);

  py::enum_<BondOrder>(m, "BondOrder")
      .value("Undefined", BondOrder::UNDEFINED)
      .value("Single", BondOrder::SINGLE)
      .value("Double", BondOrder::DOUBLE)
      .value("Triple", BondOrder::TRIPLE)
      .value("Quadruple", BondOrder::QUADRUPLE)
      .value("Aromatic", BondOrder::AROMATIC)
      .value("OneAndAHalf", BondOrder::ONEANDAHALF)
      .value("TwoAndAHalf", BondOrder::TWOANDAHALF);

  py::enum_<BondStereo>(m, "BondStereo")
      .value("Undefined", BondStereo::UNDEFINED)
      .value("None", BondStereo::NONE)
      .value("E", BondStereo::E)
      .value("Z", BondStereo::Z);

  py::enum_<ResidueType>(m, "ResidueType")
      .value("Unspecified", ResidueType::Unspecified)
      .value("AminoAcid", ResidueType::AminoAcid)
      .value("Sugar", ResidueType::Sugar)
      .value("NucleicAcid", ResidueType::NucleicAcid)
      .value("Lipid", ResidueType::Lipid)
      .value("NonSpecific", ResidueType::NonSpecific);

  // ===========================================================================
  // == Coordinates bindings ===================================================
  // ===========================================================================
  py::class_<Coordinates>(m, "Coordinates")
      .def_readonly("x", &Coordinates::x)
      .def_readonly("y", &Coordinates::y)
      .def_readonly("z", &Coordinates::z);

  // ===========================================================================
  // == Atom class bindings ====================================================
  // ===========================================================================
  py::class_<Atom>(m, "Atom")
      .def("NumBonds", &Atom::NumBonds)
      .def("NumAngles", &Atom::NumAngles)
      .def("NumDihedrals", &Atom::NumDihedrals)
      .def("NumHydrogenBonds", &Atom::NumHydrogenBonds)
      .def("NumHeavyAtomBonds", &Atom::NumHeavyAtomBonds)
      .def("HasType", &Atom::HasType)
      .def("GetElement", &Atom::GetElement)
      .def("GetFormalCharge", &Atom::GetFormalCharge)
      .def("GetPartialCharge", &Atom::GetPartialCharge)
      .def("GetTag", &Atom::GetTag)
      .def("GetChargeGroupID", &Atom::GetChargeGroupID)
      .def("GetID", &Atom::GetID)
      .def("GetResidueID", &Atom::GetResidueID)
      .def("GetResidueName", &Atom::GetResidueName)
      .def("SetResidueName", &Atom::SetResidueName)
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
      .def("SetElement", py::overload_cast<const Element &>(&Atom::SetElement))
      .def("SetElement", py::overload_cast<std::string>(&Atom::SetElement))
      .def("SetElement", py::overload_cast<int32_t>(&Atom::SetElement))
      .def("SetFormalCharge", &Atom::SetFormalCharge)
      .def("SetPartialCharge", &Atom::SetPartialCharge)
      .def("SetImplicitCount", &Atom::SetImplicitCount)
      .def("SetTag", &Atom::SetTag)
      .def("SetChargeGroupID", &Atom::SetChargeGroupID)
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
      .def("__repr__", &outstream_operator<Atom>);

  // ===========================================================================
  // == Bond class bindings ====================================================
  // ===========================================================================

  py::class_<Bond>(m, "Bond")
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
      .def("IsAmideBond", &Bond::IsAmideBond)
      .def("IsCarbonylBond", &Bond::IsCarbonylBond)
      .def("Length", &Bond::Length)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def("__bool__", &Bond::operator bool)
      .def("__str__", &outstream_operator<Bond>)
      .def("__repr__", &outstream_operator<Bond>);

  // ===========================================================================
  // == Angle class bindings ===================================================
  // ===========================================================================

  py::class_<Angle>(m, "Angle")
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
      .def("__repr__", &outstream_operator<Angle>);

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
      .def("__repr__", &outstream_operator<Dihedral>);

  // ===========================================================================
  // == Residue class bindings =================================================
  // ===========================================================================

  py::class_<Residue>(m, "Residue")
      .def("HasAtom", &Residue::HasAtom)
      .def("GetType", &Residue::GetType)
      .def("IsAminoAcid", &Residue::IsAminoAcid)
      .def("IsAlphaAminoAcid", &Residue::IsAlphaAminoAcid)
      .def("IsBetaAminoAcid", &Residue::IsBetaAminoAcid)
      .def("IsGammaAminoAcid", &Residue::IsGammaAminoAcid)
      .def("IsDeltaAminoAcid", &Residue::IsDeltaAminoAcid)
      .def("IsSugar", &Residue::IsSugar)
      .def("IsLipid", &Residue::IsLipid)
      .def("IsNucleicAcid", &Residue::IsNucleicAcid)
      .def("GetAtoms", [](const Residue &res) -> std::vector<Atom> {
        return std::vector<Atom>(res.GetAtoms().begin(), res.GetAtoms().end());
      });

  // ===========================================================================
  // == Molecule class bindings ================================================
  // ===========================================================================
  py::class_<Molecule>(m, "Molecule")
      .def(py::init<std::string>())
      .def("HasAtom", &Molecule::HasAtom)
      .def("HasBond",
           py::overload_cast<const Bond &>(&Molecule::HasBond, py::const_))
      .def("HasBond", py::overload_cast<const Atom &, const Atom &>(
                          &Molecule::HasBond, py::const_))
      .def("HasAngle",
           py::overload_cast<const Angle &>(&Molecule::HasAngle, py::const_))
      .def("HasAngle",
           py::overload_cast<const Atom &, const Atom &, const Atom &>(
               &Molecule::HasAngle))
      .def("HasDihedral", py::overload_cast<const Dihedral &>(
                              &Molecule::HasDihedral, py::const_))
      .def("HasDihedral",
           py::overload_cast<const Atom &, const Atom &, const Atom &,
                             const Atom &>(&Molecule::HasDihedral))
      .def("NumAtoms", &Molecule::NumAtoms)
      .def("NumBonds", &Molecule::NumBonds)
      .def("NumAngles", &Molecule::NumAngles)
      .def("NumDihedrals", &Molecule::NumDihedrals)
      .def("GetAtom", &Molecule::GetAtom)
      .def("GetAtomID", &Molecule::GetAtomID)
      .def("GetAtomTag", &Molecule::GetAtomTag)
      .def("GetBond",
           py::overload_cast<uint32_t>(&Molecule::GetBond, py::const_))
      .def("GetBond", py::overload_cast<const Atom &, const Atom &>(
                          &Molecule::GetBond, py::const_))
      .def("GetBondID", &Molecule::GetBondID)
      .def("GetBondTag", &Molecule::GetBondTag)
      .def("GetAngle", py::overload_cast<uint32_t>(&Molecule::GetAngle))
      .def("GetAngle",
           py::overload_cast<const Atom &, const Atom &, const Atom &>(
               &Molecule::GetAngle))
      .def("GetAngleID", &Molecule::GetAngleID)
      .def("GetAngleTag", &Molecule::GetAngleTag)
      .def("GetDihedral", py::overload_cast<uint32_t>(&Molecule::GetDihedral))
      .def("GetDihedral",
           py::overload_cast<const Atom &, const Atom &, const Atom &,
                             const Atom &>(&Molecule::GetDihedral))
      .def("GetDihedralID", &Molecule::GetDihedralID)
      .def("GetDihedralTag", &Molecule::GetDihedralTag)
      .def("GetFormula", &Molecule::GetFormula)
      .def("GetGraph", &Molecule::GetGraph)
      .def("GetCondensedGraph", &Molecule::GetCondensedGraph)
      .def("GetName", &Molecule::GetName)
      .def("GetMolecularCharge", &Molecule::GetMolecularCharge)
      .def("SetName", &Molecule::SetName)
      .def("SetMolecularCharge", &Molecule::SetMolecularCharge)

      .def("NewAtom", py::overload_cast<>(&Molecule::NewAtom))
      .def("NewAtom", py::overload_cast<const Element &>(&Molecule::NewAtom))
      .def("NewAtom",
           py::overload_cast<const Element &, double, double, double>(
               &Molecule::NewAtom))
      .def("NewBond", &Molecule::NewBond)
      .def("NewDihedral",
           py::overload_cast<const Atom &, const Atom &, const Atom &,
                             const Atom &>(&Molecule::NewDihedral))
      .def("RemoveAtom", &Molecule::RemoveAtom)
      .def("RemoveBond", py::overload_cast<const Bond &>(&Molecule::RemoveBond))
      .def("RemoveBond",
           py::overload_cast<const Atom &, const Atom &>(&Molecule::RemoveBond))
      .def("PerceiveAngles", &Molecule::PerceiveAngles)
      .def("PerceiveDihedrals", &Molecule::PerceiveDihedrals)
      .def("PerceiveResidues", &Molecule::PerceiveResidues)
      .def("OptimiseChargeGroups", &Molecule::OptimiseChargeGroups)
      .def("ReorderAtoms", &Molecule::ReorderAtoms)
      .def("UniquifyAtomNames", &Molecule::UniquifyAtomNames)
  .def("GiveAromaticBondsImpropers", &Molecule::GiveAromaticBondsImpropers)
      .def("ReserveAtoms", &Molecule::ReserveAtoms)
      .def("ReserveBonds", &Molecule::ReserveBonds)
      .def("GetAtoms", &Molecule::GetAtoms, Ref)
      .def("GetBonds", &Molecule::GetBonds, Ref)
      .def("GetAngles", &Molecule::GetAngles, Ref)
      .def("GetDihedrals", &Molecule::GetDihedrals, Ref)
      .def("GetResidueID", &Molecule::GetResidueID, Ref)
      .def("GetResidues", &Molecule::GetResidues, Ref)
      .def("GetForcefield", &Molecule::GetForcefield)
      .def("SetForcefield", &Molecule::SetForcefield)
      .def("ResetForcefield", &Molecule::ResetForcefield)
      .def("HasForcefield", &Molecule::HasForcefield)
      .def("ModificationMade", &Molecule::ModificationMade);

  // ===========================================================================
  // == Module function bindings ===============================================
  // ===========================================================================
  m.def("SaveMolecule", &SaveMolecule);
  m.def("LoadMolecule", &LoadMolecule);
  m.def("ProtonatedLysine", &special::ProtonatedLysine);
  m.def("Methionine", &special::Methionine);
  
  // Container bindings
  py::bind_vector<std::vector<Atom>>(m, "VecAtom");
  py::bind_vector<std::vector<Bond>>(m, "VecBond");
  py::bind_vector<std::vector<Angle>>(m, "VecAngle");
  py::bind_vector<std::vector<Dihedral>>(m, "VecDihedral");
}
