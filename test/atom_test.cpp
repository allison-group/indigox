#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/output_test_stream.hpp>

#include "class_test_wrappers.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>


namespace indigox::test {
  struct AtomFixture {
    Molecule mol;
    indigox::test::IXAtom a1, a2;
    Atom a3;
    Element e1, e2, enull;
    size_ num_bonds, num_angles, num_dihedrals;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    AtomFixture() : mol(new IXMolecule()), a1(Molecule()), a2(mol), a3(indigox::test::IXAtom::GetNewAtom()) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<size_> dist(1,10);
      std::array<size_, 118> elems;
      std::iota(elems.begin(), elems.end(), 1);
      std::shuffle(elems.begin(), elems.end(), gen);
      e1 = GetPeriodicTable()->GetElement(elems[0]);
      e2 = GetPeriodicTable()->GetElement(elems[1]);
      a3->SetElement(elems[2]);
      enull = GetPeriodicTable()->GetUndefined();
      num_bonds = dist(gen); bonds.reserve(num_bonds);
      for (size_ i = 0; i < num_bonds; ++i) {
        bonds.emplace_back(new indigox::IXBond());
        a2.AddBond(bonds.back());
      }
      num_angles = dist(gen); angles.reserve(num_angles);
      for (size_ i = 0; i < num_angles; ++i) {
        angles.emplace_back(new indigox::IXAngle());
        a2.AddAngle(angles.back());
      }
      num_dihedrals = dist(gen); dihedrals.reserve(num_dihedrals);
      for (size_ i = 0; i < num_dihedrals; ++i) {
        dihedrals.emplace_back(new indigox::IXDihedral());
        a2.AddDihedral(dihedrals.back());
      }
    }
  };
}

BOOST_FIXTURE_TEST_SUITE(ixatom, indigox::test::AtomFixture);

using namespace indigox;

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK(a1.GetUniqueID());
  BOOST_CHECK(a1.GetMolecule() == Molecule());
  BOOST_CHECK(a2.GetUniqueID());
  BOOST_CHECK(a1.GetUniqueID() != a2.GetUniqueID());
  BOOST_CHECK(a2.GetMolecule() == mol);
}

BOOST_AUTO_TEST_CASE(aromaticity_get_set) {
  BOOST_CHECK(!a1.GetAromaticity());  // Check init
  a1.SetAromaticity(true);
  BOOST_CHECK(a1.GetAromaticity());
  a1.SetAromaticity(false);
  BOOST_CHECK(!a1.GetAromaticity());
}

BOOST_AUTO_TEST_CASE(element_get_set) {
  // Check init to null
  BOOST_CHECK(a1.GetElement());
  BOOST_CHECK(!*a1.GetElement());
  // Check set is correct
  BOOST_CHECK_NO_THROW(a1.SetElement(e1));
  BOOST_CHECK_NO_THROW(a2.SetElement(e2));
  BOOST_CHECK(a1.GetElement() == e1);
  BOOST_CHECK(a2.GetElement() == e2);
  // Check set by name is correct
  BOOST_CHECK_NO_THROW(a1.SetElement(e2->GetName()));
  BOOST_CHECK_NO_THROW(a2.SetElement(e1->GetName()));
  BOOST_CHECK(a1.GetElement() == e2);
  BOOST_CHECK(a2.GetElement() == e1);
  // Check set by symbol is correct
  BOOST_CHECK_NO_THROW(a1.SetElement(e1->GetSymbol()));
  BOOST_CHECK_NO_THROW(a2.SetElement(e2->GetSymbol()));
  BOOST_CHECK(a1.GetElement() == e1);
  BOOST_CHECK(a2.GetElement() == e2);
  // Check set by Z is correct
  BOOST_CHECK_NO_THROW(a1.SetElement(e2->GetAtomicNumber()));
  BOOST_CHECK_NO_THROW(a2.SetElement(e1->GetAtomicNumber()));
  BOOST_CHECK(a1.GetElement() == e2);
  BOOST_CHECK(a2.GetElement() == e1);
}

BOOST_AUTO_TEST_CASE(formal_charge_get_set) {
  BOOST_CHECK(a1.GetFormalCharge() == 0);
  a2.SetFormalCharge(12);
  BOOST_CHECK(a2.GetFormalCharge() == 12);
  a1.SetFormalCharge(-14);
  BOOST_CHECK(a1.GetFormalCharge() == -14);
}

BOOST_AUTO_TEST_CASE(implicit_h_get_set) {
  BOOST_CHECK(a2.GetImplicitCount() == 0);
  a1.SetImplicitCount(3);
  BOOST_CHECK(a1.GetImplicitCount() == 3);
  a2.SetImplicitCount(0);
  BOOST_CHECK(a2.GetImplicitCount() == 0);
}

BOOST_AUTO_TEST_CASE(implicit_h_add_remove) {
  a1.SetImplicitCount(5);
  BOOST_CHECK(a1.AddImplicitHydrogen() == 6);
  BOOST_CHECK(a1.RemoveImplicitHydrogen() == 5);
  a2.SetImplicitCount(1);
  BOOST_CHECK(a2.RemoveImplicitHydrogen() == 0);
  BOOST_CHECK(a2.RemoveImplicitHydrogen() == 0);
  BOOST_CHECK(a2.AddImplicitHydrogen() == 1);
}

BOOST_AUTO_TEST_CASE(molecule_get) {
  BOOST_CHECK(a1.GetMolecule() == Molecule());
  BOOST_CHECK(a2.GetMolecule() == mol);
  // Should return empty Molecule if assigned molecule is deleted
  mol.reset();
  BOOST_CHECK(!a2.GetMolecule());
  BOOST_CHECK(a2.GetMolecule() == Molecule());
}

BOOST_AUTO_TEST_CASE(name_get_set) {
  BOOST_CHECK(a1.GetName() == "");
  a1.SetName("TEST");
  BOOST_CHECK(a1.GetName() == "TEST");
  a2.SetName("MorE");
  BOOST_CHECK(a2.GetName() == "MorE");
}

BOOST_AUTO_TEST_CASE(partial_charge_get_set) {
  BOOST_CHECK_CLOSE(a1.GetPartialCharge(), static_cast<float_>(0.0), 0.0001);
  a1.SetPartialCharge(0.7785);
  BOOST_CHECK_CLOSE(a1.GetPartialCharge(), static_cast<float_>(0.7785), 0.0001);
  a2.SetPartialCharge(-7.359);
  BOOST_CHECK_CLOSE(a2.GetPartialCharge(), static_cast<float_>(-7.359), 0.0001);
}

BOOST_AUTO_TEST_CASE(position_get_set) {
  BOOST_CHECK_CLOSE(a1.GetX(), static_cast<float_>(0.0), 0.0001);
  BOOST_CHECK_CLOSE(a1.GetY(), static_cast<float_>(0.0), 0.0001);
  BOOST_CHECK_CLOSE(a1.GetZ(), static_cast<float_>(0.0), 0.0001);
  a1.SetPosition(1.1, -1.5, 12.89);
  Vec3 pos = a1.GetVector();
  BOOST_CHECK_CLOSE(a1.GetX(), static_cast<float_>(1.1), 0.0001);
  BOOST_CHECK_CLOSE(a1.GetY(), static_cast<float_>(-1.5), 0.0001);
  BOOST_CHECK_CLOSE(a1.GetZ(), static_cast<float_>(12.89), 0.0001);
  BOOST_CHECK_CLOSE(pos.x, static_cast<float_>(1.1), 0.0001);
  BOOST_CHECK_CLOSE(pos.y, static_cast<float_>(-1.5), 0.0001);
  BOOST_CHECK_CLOSE(pos.z, static_cast<float_>(12.89), 0.0001);
  
  a2.SetX(-72.00001);
  a2.SetY(-0.000009);
  a2.SetZ(1.12345);
  BOOST_CHECK_CLOSE(a2.GetX(), static_cast<float_>(-72.00001), 0.0001);
  BOOST_CHECK_CLOSE(a2.GetY(), static_cast<float_>(-0.000009), 0.0001);
  BOOST_CHECK_CLOSE(a2.GetZ(), static_cast<float_>(1.12345), 0.0001);
  pos = a2.GetVector();
  BOOST_CHECK_CLOSE(pos.x, static_cast<float_>(-72.00001), 0.0001);
  BOOST_CHECK_CLOSE(pos.y, static_cast<float_>(-0.000009), 0.0001);
  BOOST_CHECK_CLOSE(pos.z, static_cast<float_>(1.12345), 0.0001);
}

BOOST_AUTO_TEST_CASE(stereochemistry_get_set) {
  BOOST_CHECK(a1.GetStereochemistry() == AtomStereo::UNDEFINED);
  a1.SetStereochemistry(AtomStereo::ACHIRAL);
  BOOST_CHECK(a1.GetStereochemistry() == AtomStereo::ACHIRAL);
  a2.SetStereochemistry(AtomStereo::S);
  BOOST_CHECK(a2.GetStereochemistry() == AtomStereo::S);
  a1.SetStereochemistry(AtomStereo::R);
  BOOST_CHECK(a1.GetStereochemistry() == AtomStereo::R);
}

BOOST_AUTO_TEST_CASE(tag_get_set) {
  BOOST_CHECK(a1.GetTag() == 0);
  a1.SetTag(12);
  BOOST_CHECK(a1.GetTag() == 12);
  a1.SetTag(0);
  BOOST_CHECK(a1.GetTag() == 0);
}

BOOST_AUTO_TEST_CASE(printing_methods) {
  // Test ostream operator
  boost::test_tools::output_test_stream actual_output;
  std::stringstream expected_out;
  expected_out << "Atom(" << a3->GetUniqueID() << ", " << a3->GetElement()->GetSymbol() << ")";
  actual_output << a3;
  BOOST_CHECK(actual_output.is_equal(expected_out.str()));
  actual_output << Atom();
  BOOST_CHECK(actual_output.is_empty());
  
  // Test ToString method
  BOOST_CHECK(a1.ToString() == "Atom(, XX)");
  a2.SetName("Longnametest w i t h s p a c e s");
  a2.SetElement(34);
  BOOST_CHECK(a2.ToString() == "Atom(Longnametest w i t h s p a c e s, Se)");
  a1.SetName("test");
  a1.SetElement(118);
  BOOST_CHECK(a1.ToString() == "Atom(test, Og)");
}

BOOST_AUTO_TEST_CASE(angles_add_remove) {
  BOOST_CHECK(a1.NumAngles() == 0);
  BOOST_CHECK(a2.NumAngles() == num_angles);
  // Check iterators
  auto angle_it = a2.GetAngleIters();
  BOOST_CHECK(std::distance(angle_it.first, angle_it.second) == num_angles);
  std::vector<Angle> stored_angles;
  for (; angle_it.first != angle_it.second; ++angle_it.first)
    stored_angles.emplace_back(angle_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(angles.begin(), angles.end(),
                                stored_angles.begin(), stored_angles.end());
  // Remove the last angle
  a2.RemoveAngle(angles.back());
  BOOST_CHECK(a2.NumAngles() == num_angles - 1);
  angle_it = a2.GetAngleIters();
  BOOST_CHECK(std::distance(angle_it.first, angle_it.second) == num_angles - 1);
  stored_angles.clear(); angles.resize(num_angles - 1);
  for (; angle_it.first != angle_it.second; ++angle_it.first)
    stored_angles.emplace_back(angle_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(angles.begin(), angles.end(),
                                stored_angles.begin(), stored_angles.end());
}

BOOST_AUTO_TEST_CASE(bonds_add_remove) {
  BOOST_CHECK(a1.NumBonds() == 0);
  BOOST_CHECK(a2.NumBonds() == num_bonds);
  // Check iterators
  auto bond_it = a2.GetBondIters();
  BOOST_CHECK(std::distance(bond_it.first, bond_it.second) == num_bonds);
  std::vector<Bond> stored_bonds;
  for (; bond_it.first != bond_it.second; ++bond_it.first)
    stored_bonds.emplace_back(bond_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(bonds.begin(), bonds.end(),
                                stored_bonds.begin(), stored_bonds.end());
  // Remove the last bond
  a2.RemoveBond(bonds.back());
  BOOST_CHECK(a2.NumBonds() == num_bonds - 1);
  bond_it = a2.GetBondIters();
  BOOST_CHECK(std::distance(bond_it.first, bond_it.second) == num_bonds - 1);
  stored_bonds.clear(); bonds.resize(num_bonds - 1);
  for (; bond_it.first != bond_it.second; ++bond_it.first)
    stored_bonds.emplace_back(bond_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(bonds.begin(), bonds.end(),
                                stored_bonds.begin(), stored_bonds.end());
}

BOOST_AUTO_TEST_CASE(dihedrals_add_remove) {
  BOOST_CHECK(a1.NumDihedrals() == 0);
  BOOST_CHECK(a2.NumDihedrals() == num_dihedrals);
  // Check iterators
  auto dihedral_it = a2.GetDihedralIters();
  BOOST_CHECK(std::distance(dihedral_it.first, dihedral_it.second) == num_dihedrals);
  std::vector<Dihedral> stored_dihedrals;
  for (; dihedral_it.first != dihedral_it.second; ++dihedral_it.first)
    stored_dihedrals.emplace_back(dihedral_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(dihedrals.begin(), dihedrals.end(),
                                stored_dihedrals.begin(), stored_dihedrals.end());
  // Remove the last dihedral
  a2.RemoveDihedral(dihedrals.back());
  BOOST_CHECK(a2.NumDihedrals() == num_dihedrals - 1);
  dihedral_it = a2.GetDihedralIters();
  BOOST_CHECK(std::distance(dihedral_it.first, dihedral_it.second) == num_dihedrals - 1);
  stored_dihedrals.clear(); dihedrals.resize(num_dihedrals - 1);
  for (; dihedral_it.first != dihedral_it.second; ++dihedral_it.first)
    stored_dihedrals.emplace_back(dihedral_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(dihedrals.begin(), dihedrals.end(),
                                stored_dihedrals.begin(), stored_dihedrals.end());
}

BOOST_AUTO_TEST_CASE(cleaning_methods) {
  // Clear everything
  a2.Clear();
  BOOST_CHECK(a2.GetAromaticity() == false);
  BOOST_CHECK(a2.GetElement());
  BOOST_CHECK(!*a2.GetElement());
  BOOST_CHECK(a2.GetFormalCharge() == 0);
  BOOST_CHECK(a2.GetImplicitCount() == 0);
  BOOST_CHECK(a2.GetMolecule() == Molecule());
  BOOST_CHECK(a2.GetName() == "");
  BOOST_CHECK_CLOSE(a2.GetPartialCharge(), static_cast<float_>(0.0), 0.0001);
  BOOST_CHECK(a2.GetStereochemistry() == AtomStereo::UNDEFINED);
  BOOST_CHECK(a2.GetTag() == 0);
  BOOST_CHECK_CLOSE(a2.GetX(), static_cast<float_>(0.0), 0.0001);
  BOOST_CHECK_CLOSE(a2.GetY(), static_cast<float_>(0.0), 0.0001);
  BOOST_CHECK_CLOSE(a2.GetZ(), static_cast<float_>(0.0), 0.0001);
  BOOST_CHECK(a2.NumAngles() == 0);
  BOOST_CHECK(a2.NumBonds() == 0);
  BOOST_CHECK(a2.NumDihedrals() == 0);
}

BOOST_AUTO_TEST_SUITE_END();
