#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/helpers.hpp>

#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace indigox {
  // Dummyclasses
  class IXAngle {};
  class IXDihedral {};
  class IXMolecule {};
}

BOOST_AUTO_TEST_SUITE(ixatom);

using namespace indigox;

BOOST_AUTO_TEST_CASE(constructor) {
  IXAtom a1 = IXAtom();
  BOOST_TEST(a1.GetUniqueID());
  BOOST_TEST(!a1.GetMolecule());
  
  Molecule mol = Molecule(new IXMolecule());
  IXAtom a2 = IXAtom(mol);
  BOOST_TEST(a1.GetUniqueID());
  BOOST_TEST((a2.GetMolecule() && a2.GetMolecule() == mol));
}

BOOST_AUTO_TEST_CASE(aromaticity_get_set) {
  IXAtom a = IXAtom();
  BOOST_TEST(!a.GetAromaticity());  // Check init
  a.SetAromaticity(true);
  BOOST_TEST(a.GetAromaticity());
  a.SetAromaticity(false);
  BOOST_TEST(!a.GetAromaticity());
}

BOOST_AUTO_TEST_CASE(element_get_set) {
  PeriodicTable pt = GetPeriodicTable();
  Element carbon = pt->GetElement(6), rubidium = pt->GetElement(37);
  IXAtom a = IXAtom();
  BOOST_TEST((a.GetElement() && !*a.GetElement()));  // Check init
  BOOST_TEST(a.GetElement()->GetSymbol() == "XX");
  BOOST_CHECK_NO_THROW(a.SetElement(carbon));
  BOOST_TEST((a.GetElement() == carbon && a.GetElement() != rubidium));
  BOOST_CHECK_NO_THROW(a.SetElement(rubidium));
  BOOST_TEST((a.GetElement() == rubidium && a.GetElement() != carbon));
  BOOST_CHECK_NO_THROW(a.SetElement("C"));
  BOOST_TEST((a.GetElement() == carbon && a.GetElement() != rubidium));
  BOOST_CHECK_NO_THROW(a.SetElement("Rb"));
  BOOST_TEST((a.GetElement() == rubidium && a.GetElement() != carbon));
  BOOST_CHECK_NO_THROW(a.SetElement("carBON"));
  BOOST_TEST((a.GetElement() == carbon && a.GetElement() != rubidium));
  BOOST_CHECK_NO_THROW(a.SetElement("Rubidium"));
  BOOST_TEST((a.GetElement() == rubidium && a.GetElement() != carbon));
  BOOST_CHECK_NO_THROW(a.SetElement(6));
  BOOST_TEST((a.GetElement() == carbon && a.GetElement() != rubidium));
  BOOST_CHECK_NO_THROW(a.SetElement(37));
  BOOST_TEST((a.GetElement() == rubidium && a.GetElement() != carbon));
  // Calls that should throw errors
  BOOST_CHECK_THROW(a.SetElement(0), std::invalid_argument);
  BOOST_CHECK_THROW(a.SetElement(-1), std::invalid_argument);
  BOOST_CHECK_THROW(a.SetElement("c"), std::invalid_argument);
  BOOST_CHECK_THROW(a.SetElement("notARealElementium"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(formal_charge_get_set) {
  IXAtom a = IXAtom();
  BOOST_TEST(a.GetFormalCharge() == 0);
  a.SetFormalCharge(12);
  BOOST_TEST(a.GetFormalCharge() == 12);
  a.SetFormalCharge(-14);
  BOOST_TEST(a.GetFormalCharge() == -14);
}

BOOST_AUTO_TEST_CASE(implicit_h_get_set) {
  IXAtom a = IXAtom();
  BOOST_TEST(a.GetImplicitCount() == 0);
  a.SetImplicitCount(3);
  BOOST_TEST(a.GetImplicitCount() == 3);
  a.SetImplicitCount(0);
  BOOST_TEST(a.GetImplicitCount() == 0);
}

BOOST_AUTO_TEST_CASE(implicit_h_add_remove) {
  IXAtom a = IXAtom();
  a.SetImplicitCount(5);
  BOOST_TEST(a.AddImplicitHydrogen() == 6);
  BOOST_TEST(a.RemoveImplicitHydrogen() == 5);
  a.SetImplicitCount(1);
  BOOST_TEST(a.RemoveImplicitHydrogen() == 0);
  BOOST_TEST(a.RemoveImplicitHydrogen() == 0);
  BOOST_TEST(a.AddImplicitHydrogen() == 1);
}

BOOST_AUTO_TEST_CASE(molecule_get_set) {
  Molecule mol1 = Molecule(new IXMolecule());
  Molecule mol2 = Molecule(new IXMolecule());
  IXAtom a = IXAtom(mol1);
  BOOST_TEST(a.GetMolecule() == mol1);
  a.SetMolecule(mol2);
  BOOST_TEST((a.GetMolecule() == mol2 && a.GetMolecule() != mol1));
  // Should return empty Molecule if assigned molecule is deleted
  mol2.reset();
  BOOST_TEST(!a.GetMolecule());
  a.SetMolecule(mol1);
  BOOST_TEST(a.GetMolecule() == mol1);
}

BOOST_AUTO_TEST_CASE(name_get_set) {
  IXAtom a = IXAtom();
  BOOST_TEST(a.GetName() == "");
  a.SetName("TEST");
  BOOST_TEST(a.GetName() == "TEST");
  a.SetName("MorE");
  BOOST_TEST(a.GetName() == "MorE");
}

BOOST_AUTO_TEST_CASE(partial_charge_get_set) {
  IXAtom a = IXAtom();
  BOOST_TEST(a.GetPartialCharge() == 0.0);
  a.SetPartialCharge(0.7785);
  BOOST_TEST(a.GetPartialCharge() == 0.7785,
             boost::test_tools::tolerance(0.00001));
  a.SetPartialCharge(-7.359);
  BOOST_TEST(a.GetPartialCharge() == -7.359,
             boost::test_tools::tolerance(0.00001));
}

BOOST_AUTO_TEST_CASE(position_get_set) {
  IXAtom a = IXAtom();
  Vec3 pos = a.GetVector();
  BOOST_TEST((pos.x == 0.0 && pos.y == 0.0 && pos.z == 0.0));
  BOOST_TEST(a.GetX() == 0.0);
  BOOST_TEST(a.GetY() == 0.0);
  BOOST_TEST(a.GetZ() == 0.0);
  a.SetPosition(1.1, -1.5, 12.89);
  pos = a.GetVector();
  BOOST_TEST((pos.x == 1.1 && pos.y == -1.5 && pos.z == 12.89),
             boost::test_tools::tolerance(0.00001));
  BOOST_TEST(a.GetX() == 1.1, boost::test_tools::tolerance(0.00001));
  BOOST_TEST(a.GetY() == -1.5, boost::test_tools::tolerance(0.00001));
  BOOST_TEST(a.GetZ() == 12.89, boost::test_tools::tolerance(0.00001));
  a.SetX(-72.00001);
  BOOST_TEST(a.GetX() == -72.00001, boost::test_tools::tolerance(0.00001));
  a.SetY(-0.000009);
  BOOST_TEST(a.GetY() == -0.000009, boost::test_tools::tolerance(0.00001));
  a.SetZ(1.12345);
  BOOST_TEST(a.GetZ() == 1.12345, boost::test_tools::tolerance(0.00001));
  pos = a.GetVector();
  BOOST_TEST((pos.x == -72.00001 && pos.y == -0.000009 && pos.z == 1.12345),
             boost::test_tools::tolerance(0.0000001));
  
}

BOOST_AUTO_TEST_CASE(stereochemistry_get_set) {
  IXAtom a = IXAtom();
  BOOST_CHECK(a.GetStereochemistry() == IXAtom::Stereo::UNDEFINED);
  a.SetStereochemistry(IXAtom::Stereo::ACHIRAL);
  BOOST_CHECK(a.GetStereochemistry() == IXAtom::Stereo::ACHIRAL);
  a.SetStereochemistry(IXAtom::Stereo::S);
  BOOST_CHECK(a.GetStereochemistry() == IXAtom::Stereo::S);
  a.SetStereochemistry(IXAtom::Stereo::R);
  BOOST_CHECK(a.GetStereochemistry() == IXAtom::Stereo::R);
}

BOOST_AUTO_TEST_CASE(tag_get_set) {
  IXAtom a = IXAtom();
  BOOST_TEST(a.GetTag() == 0);
  a.SetTag(12);
  BOOST_TEST(a.GetTag() == 12);
  a.SetTag(0);
  BOOST_TEST(a.GetTag() == 0);
}

BOOST_AUTO_TEST_CASE(printing_methods) {
  // Test ToString
  IXAtom a = IXAtom();
  BOOST_TEST(a.ToString() == "Atom(, XX)");
  a.SetName("Longnametest w i t h s p a c e s");
  a.SetElement(34);
  BOOST_TEST(a.ToString() == "Atom(Longnametest w i t h s p a c e s, Se)");
  a.SetName("test");
  a.SetElement(118);
  BOOST_TEST(a.ToString() == "Atom(test, Og)");
}

BOOST_AUTO_TEST_CASE(angles_add_remove) {
  IXAtom a = IXAtom();
  Angle ang1 = Angle(new IXAngle());
  Angle ang2 = Angle(new IXAngle());
  Angle ang3 = Angle(new IXAngle());
  Angle ang4 = Angle(new IXAngle());
  Angle ang5 = Angle(new IXAngle());
  BOOST_TEST(a.NumAngles() == 0);
  // Add a single angle
  BOOST_TEST(a.AddAngle(ang1));
  BOOST_TEST(a.NumAngles() == 1);
  // Add another angle
  BOOST_TEST(a.AddAngle(ang2));
  BOOST_TEST(a.NumAngles() == 2);
  // Check iterators
  auto angles = a.GetAngleIters();
  BOOST_CHECK(std::distance(angles.first, angles.second) == 2);
  // Add an existing angle
  BOOST_TEST(!a.AddAngle(ang1));
  BOOST_TEST(a.NumAngles() == 2);
  // Remove an existing angle
  BOOST_TEST(a.RemoveAngle(ang1));
  BOOST_TEST(a.NumAngles() == 1);
  // Attempt to remove angles that no longer exist or never existed
  BOOST_TEST((!a.RemoveAngle(ang1) && !a.RemoveAngle(ang3)));
  // Add in remaining angles
  BOOST_TEST((a.AddAngle(ang3) && a.AddAngle(ang4) && a.AddAngle(ang5)));
  BOOST_TEST(a.NumAngles() == 4);
  // Adding/removing empty angle should fail
  BOOST_TEST((!a.AddAngle(Angle()) && !a.RemoveAngle(Angle())));
  // Expired angles shouldn't be counted
  ang5.reset(); ang2.reset();
  BOOST_TEST(a.NumAngles() == 2);
}

BOOST_AUTO_TEST_CASE(bonds_add_remove) {
  IXAtom a = IXAtom();
  Bond bnd1 = Bond(new IXBond());
  Bond bnd2 = Bond(new IXBond());
  Bond bnd3 = Bond(new IXBond());
  Bond bnd4 = Bond(new IXBond());
  Bond bnd5 = Bond(new IXBond());
  BOOST_TEST(a.NumBonds() == 0);
  // Add a single bond
  BOOST_TEST(a.AddBond(bnd1));
  BOOST_TEST(a.NumBonds() == 1);
  // Add another bond
  BOOST_TEST(a.AddBond(bnd2));
  BOOST_TEST(a.NumBonds() == 2);
  // Check iterators
  auto bonds = a.GetBondIters();
  BOOST_CHECK(std::distance(bonds.first, bonds.second) == 2);
  // Add an existing bond
  BOOST_TEST(!a.AddBond(bnd1));
  BOOST_TEST(a.NumBonds() == 2);
  // Remove an existing bond
  BOOST_TEST(a.RemoveBond(bnd1));
  BOOST_TEST(a.NumBonds() == 1);
  // Attempt to remove bonds that no longer exist or never existed
  BOOST_TEST((!a.RemoveBond(bnd1) && !a.RemoveBond(bnd3)));
  // Add in remaining bonds
  BOOST_TEST((a.AddBond(bnd3) && a.AddBond(bnd4) && a.AddBond(bnd5)));
  BOOST_TEST(a.NumBonds() == 4);
  // Adding/removing empty bond should fail
  BOOST_TEST((!a.AddBond(Bond()) && !a.RemoveBond(Bond())));
  // Expired bonds shouldn't be counted
  bnd5.reset(); bnd2.reset();
  BOOST_TEST(a.NumBonds() == 2);
}

BOOST_AUTO_TEST_CASE(dihedrals_add_remove) {
  IXAtom a = IXAtom();
  Dihedral dhd1 = Dihedral(new IXDihedral());
  Dihedral dhd2 = Dihedral(new IXDihedral());
  Dihedral dhd3 = Dihedral(new IXDihedral());
  Dihedral dhd4 = Dihedral(new IXDihedral());
  Dihedral dhd5 = Dihedral(new IXDihedral());
  BOOST_TEST(a.NumDihedrals() == 0);
  // Add a single dihedral
  BOOST_TEST(a.AddDihedral(dhd1));
  BOOST_TEST(a.NumDihedrals() == 1);
  // Add another dihedral
  BOOST_TEST(a.AddDihedral(dhd2));
  BOOST_TEST(a.NumDihedrals() == 2);
  // Check iterators
  auto dihedrals = a.GetDihedralIters();
  BOOST_CHECK(std::distance(dihedrals.first, dihedrals.second) == 2);
  // Add an existing dihedral
  BOOST_TEST(!a.AddDihedral(dhd1));
  BOOST_TEST(a.NumDihedrals() == 2);
  // Remove an existing dihedral
  BOOST_TEST(a.RemoveDihedral(dhd1));
  BOOST_TEST(a.NumDihedrals() == 1);
  // Attempt to remove dihedrals that no longer exist or never existed
  BOOST_TEST((!a.RemoveDihedral(dhd1) && !a.RemoveDihedral(dhd3)));
  // Add in remaining dihedrals
  BOOST_TEST((a.AddDihedral(dhd3) && a.AddDihedral(dhd4) && a.AddDihedral(dhd5)));
  BOOST_TEST(a.NumDihedrals() == 4);
  // Adding/removing empty dihedral should fail
  BOOST_TEST((!a.AddDihedral(Dihedral()) && !a.RemoveDihedral(Dihedral())));
  // Expired dihedrals shouldn't be counted
  dhd5.reset(); dhd2.reset();
  BOOST_TEST(a.NumDihedrals() == 2);
}

BOOST_AUTO_TEST_CASE(cleaning_methods) {
  // Build a resonably populated atom
  Molecule mol = Molecule(new IXMolecule());
  IXAtom a(mol);
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;
  for (int i = 0; i < 5; ++i) {
    bonds.emplace_back(new IXBond());
    a.AddBond(bonds.back());
    angles.emplace_back(new IXAngle());
    a.AddAngle(angles.back());
    dihedrals.emplace_back(new IXDihedral());
    a.AddDihedral(dihedrals.back());
  }
  a.SetElement("F");
  a.SetFormalCharge(12);
  a.SetPartialCharge(1.34);
  a.SetTag(4);
  a.SetImplicitCount(3);
  a.SetName("tester");
  a.SetPosition(3.0, 3.0, 3.0);
  a.SetStereochemistry(IXAtom::Stereo::R);
  a.SetAromaticity(true);
  // Check everything nicely set
  BOOST_TEST((a.GetAromaticity() && a.GetElement() == "F" && a.GetFormalCharge() == 12
              && a.GetImplicitCount() == 3 && a.GetMolecule() == mol && a.GetName() == "tester"
              && a.GetPartialCharge() == 1.34 && a.GetStereochemistry() == IXAtom::Stereo::R
              && a.GetTag() == 4 && a.GetX() == 3.0 && a.GetY() == 3.0 && a.GetZ() == 3.0
              && a.NumAngles() == 5 && a.NumBonds() == 5 && a.NumDihedrals() == 5),
             boost::test_tools::tolerance(0.0000001));
  // Invalidate some bonds, angles, dihedrals, then cleanup
  bonds[0].reset();
  angles[0].reset(); angles[1].reset();
  dihedrals[0].reset(); dihedrals[1].reset(); dihedrals[2].reset();
  IXAtom::AtomBondIter b_begin, b_end;
  IXAtom::AtomAngleIter a_begin, a_end;
  IXAtom::AtomDihedralIter d_begin, d_end;
  std::tie(b_begin, b_end) = a.GetBondIters();
  std::tie(a_begin, a_end) = a.GetAngleIters();
  std::tie(d_begin, d_end) = a.GetDihedralIters();
  BOOST_CHECK(std::distance(b_begin, b_end) == 5);
  BOOST_CHECK(std::distance(a_begin, a_end) == 5);
  BOOST_CHECK(std::distance(d_begin, d_end) == 5);
  a.Cleanup();
  std::tie(b_begin, b_end) = a.GetBondIters();
  std::tie(a_begin, a_end) = a.GetAngleIters();
  std::tie(d_begin, d_end) = a.GetDihedralIters();
  BOOST_CHECK(std::distance(b_begin, b_end) == 4);
  BOOST_CHECK(std::distance(a_begin, a_end) == 3);
  BOOST_CHECK(std::distance(d_begin, d_end) == 2);
  
  // Clear everything
  a.Clear();
  BOOST_CHECK(!a.GetAromaticity());
  BOOST_CHECK(!*a.GetElement());
  BOOST_CHECK(a.GetFormalCharge() == 0);
  BOOST_CHECK(a.GetImplicitCount() == 0);
  BOOST_CHECK(!a.GetMolecule());
  BOOST_CHECK(a.GetName() == "");
  BOOST_CHECK(a.GetPartialCharge() == 0.0);
  BOOST_CHECK(a.GetStereochemistry() == IXAtom::Stereo::UNDEFINED);
  BOOST_CHECK(a.GetTag() == 0);
  BOOST_CHECK((a.GetX() == 0.0 && a.GetY() == 0.0 && a.GetZ() == 0.0));
  BOOST_CHECK((a.NumAngles() == 0 && a.NumBonds() == 0 && a.NumDihedrals() == 0));
}

BOOST_AUTO_TEST_SUITE_END();
