#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

#include "class_test_wrappers.hpp"

#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace indigox::test {
  struct BondFixture {
    Molecule mol;
    Atom a1, a2, a3, a4;
    indigox::test::IXBond b1, b2;
    size_ num_angles, num_dihedrals;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    
    BondFixture() : mol(new IXMolecule()), a1(indigox::test::IXAtom::GetNewAtom()),
    a2(indigox::test::IXAtom::GetNewAtom(mol)), a3(indigox::test::IXAtom::GetNewAtom()),
    a4(indigox::test::IXAtom::GetNewAtom(mol)), b1(a1,a3,Molecule()), b2(a2,a4,mol) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<size_> dist(3,10);
      num_angles = dist(gen); num_dihedrals = dist(gen);
      for (size_ i = 0; i < num_angles; ++i) {
        angles.emplace_back(new indigox::IXAngle());
        b2.AddAngle(angles.back());
      }
      for (size_ i = 0; i < num_dihedrals; ++i) {
        dihedrals.emplace_back(new indigox::IXDihedral());
        b2.AddDihedral(dihedrals.back());
      }
    }
    
  };
}

using namespace indigox;

BOOST_FIXTURE_TEST_SUITE(ixbond, indigox::test::BondFixture);

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK(b1.GetUniqueID());
  BOOST_CHECK(b2.GetUniqueID());
  BOOST_CHECK(b1.GetUniqueID() != b2.GetUniqueID());
}

BOOST_AUTO_TEST_CASE(aromaticity_get_set) {
  BOOST_CHECK(!b1.GetAromaticity());
  b1.SetAromaticity(true);
  BOOST_CHECK(b1.GetAromaticity());
  b1.SetAromaticity(false);
  BOOST_CHECK(!b1.GetAromaticity());
}

BOOST_AUTO_TEST_CASE(molecule_get) {
  BOOST_CHECK(b1.GetMolecule() == Molecule());
  BOOST_CHECK(b2.GetMolecule() == mol);
  mol.reset();
  BOOST_CHECK(!b2.GetMolecule());
  BOOST_CHECK(b2.GetMolecule() == Molecule());
}

BOOST_AUTO_TEST_CASE(order_get_set) {
  BOOST_CHECK(b1.GetOrder() == BondOrder::UNDEFINED);
  b1.SetOrder(BondOrder::TRIPLE);
  BOOST_CHECK(b1.GetOrder() == BondOrder::TRIPLE);
  b2.SetOrder(BondOrder::TWOANDAHALF);
  BOOST_CHECK(b2.GetOrder() == BondOrder::TWOANDAHALF);
}

BOOST_AUTO_TEST_CASE(atoms_get_swap) {
  BOOST_CHECK(b1.NumAtoms() == 2);
  BOOST_CHECK(b1.GetSourceAtom() == a1);
  BOOST_CHECK(b1.GetTargetAtom() == a3);
  BOOST_CHECK(b2.GetSourceAtom() == a2);
  BOOST_CHECK(b2.GetTargetAtom() == a4);
  a1.reset();
  BOOST_CHECK(b1.GetSourceAtom() == Atom());
  BOOST_CHECK(b1.NumAtoms() == 2);
  
  b2.SwapSourceTarget();
  BOOST_CHECK(b2.GetSourceAtom() == a4);
  BOOST_CHECK(b2.GetTargetAtom() == a2);
  
  // Check iterators
  auto atm_it = b2.GetAtomIters();
  BOOST_CHECK(std::distance(atm_it.first, atm_it.second) == 2);
  std::vector<Atom> expect_atm = {a4,a2};
  std::vector<Atom> obtain_atm;
  for (; atm_it.first != atm_it.second; ++atm_it.first)
    obtain_atm.emplace_back(atm_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(expect_atm.begin(), expect_atm.end(),
                                obtain_atm.begin(), obtain_atm.end());
}

BOOST_AUTO_TEST_CASE(stereochemistry_get_set) {
  BOOST_CHECK(b1.GetStereochemistry() == BondStereo::UNDEFINED);
  b1.SetStereochemistry(BondStereo::NONE);
  BOOST_CHECK(b1.GetStereochemistry() == BondStereo::NONE);
  b2.SetStereochemistry(BondStereo::Z);
  BOOST_CHECK(b2.GetStereochemistry() == BondStereo::Z);
}

BOOST_AUTO_TEST_CASE(tag_get_set) {
  BOOST_CHECK(b1.GetTag() == 0);
  b1.SetTag(12);
  BOOST_CHECK(b1.GetTag() == 12);
  b1.SetTag(123);
  BOOST_CHECK(b1.GetTag() == 123);
}

BOOST_AUTO_TEST_CASE(angle_add_remove) {
  BOOST_CHECK(b1.NumAngles() == 0);
  BOOST_CHECK(b2.NumAngles() == num_angles);
  // Check iterators
  auto angle_it = b2.GetAngleIters();
  BOOST_CHECK(std::distance(angle_it.first, angle_it.second) == num_angles);
  std::vector<Angle> stored_angles;
  for (; angle_it.first != angle_it.second; ++angle_it.first)
    stored_angles.emplace_back(angle_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(angles.begin(), angles.end(),
                                stored_angles.begin(), stored_angles.end());
  // Remove the last angle
  b2.RemoveAngle(angles.back());
  BOOST_CHECK(b2.NumAngles() == num_angles - 1);
  angle_it = b2.GetAngleIters();
  BOOST_CHECK(std::distance(angle_it.first, angle_it.second) == num_angles - 1);
  stored_angles.clear(); angles.resize(num_angles - 1);
  for (; angle_it.first != angle_it.second; ++angle_it.first)
    stored_angles.emplace_back(angle_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(angles.begin(), angles.end(),
                                stored_angles.begin(), stored_angles.end());
}

BOOST_AUTO_TEST_CASE(dihedrals_add_remove) {
  BOOST_CHECK(b1.NumDihedrals() == 0);
  BOOST_CHECK(b2.NumDihedrals() == num_dihedrals);
  // Check iterators
  auto dihedral_it = b2.GetDihedralIters();
  BOOST_CHECK(std::distance(dihedral_it.first, dihedral_it.second) == num_dihedrals);
  std::vector<Dihedral> stored_dihedrals;
  for (; dihedral_it.first != dihedral_it.second; ++dihedral_it.first)
    stored_dihedrals.emplace_back(dihedral_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(dihedrals.begin(), dihedrals.end(),
                                stored_dihedrals.begin(), stored_dihedrals.end());
  // Remove the last dihedral
  b2.RemoveDihedral(dihedrals.back());
  BOOST_CHECK(b2.NumDihedrals() == num_dihedrals - 1);
  dihedral_it = b2.GetDihedralIters();
  BOOST_CHECK(std::distance(dihedral_it.first, dihedral_it.second) == num_dihedrals - 1);
  stored_dihedrals.clear(); dihedrals.resize(num_dihedrals - 1);
  for (; dihedral_it.first != dihedral_it.second; ++dihedral_it.first)
    stored_dihedrals.emplace_back(dihedral_it.first->lock());
  BOOST_CHECK_EQUAL_COLLECTIONS(dihedrals.begin(), dihedrals.end(),
                                stored_dihedrals.begin(), stored_dihedrals.end());
}

BOOST_AUTO_TEST_CASE(printing_methods) {
  // Test ostream operator
  boost::test_tools::output_test_stream actual_output;
  actual_output << indigox::test::IXBond::GetNewBond();
  BOOST_CHECK(actual_output.is_equal("Bond(, )"));
  
  // Test ToString
  BOOST_CHECK(b1.ToString() == "Bond(Atom(, XX), Atom(, XX))");
  a1->SetName("Romeo"); a1->SetElement("U");
  a3->SetName("Foxtrot"); a3->SetElement("F");
  b1.SwapSourceTarget();
  BOOST_CHECK(b1.ToString() == "Bond(Atom(Foxtrot, F), Atom(Romeo, U))");
}

BOOST_AUTO_TEST_CASE(cleaning_methods) {
  b2.Clear();
  BOOST_CHECK(b2.GetAromaticity() == false);
  BOOST_CHECK(b2.GetMolecule() == Molecule());
  BOOST_CHECK(b2.GetOrder() == BondOrder::UNDEFINED);
  BOOST_CHECK(b2.GetSourceAtom() == Atom());
  BOOST_CHECK(b2.GetStereochemistry() == BondStereo::UNDEFINED);
  BOOST_CHECK(b2.GetTag() == 0);
  BOOST_CHECK(b2.GetTargetAtom() == Atom());
  BOOST_CHECK(b2.NumAngles() == 0);
  BOOST_CHECK(b2.NumAtoms() == 2);
  BOOST_CHECK(b2.NumDihedrals() == 0);
}


BOOST_AUTO_TEST_SUITE_END();
