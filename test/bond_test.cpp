#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/utils/numerics.hpp>

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

using namespace indigox;

BOOST_AUTO_TEST_SUITE(ixbond);

BOOST_AUTO_TEST_CASE(constructor) {
  Atom a1 = Atom(new IXAtom());
  Atom a2 = Atom(new IXAtom());
  // Check throwing conditions
  BOOST_CHECK_NO_THROW(IXBond(a1,a2));
  BOOST_CHECK_THROW(IXBond(a1,a1), std::logic_error);
  
  IXBond b1 = IXBond();
  IXBond b2 = IXBond(a1, a2);
  // Check construction conditions
  BOOST_CHECK(b1.GetUniqueID());
  BOOST_CHECK(b2.GetUniqueID());
  BOOST_TEST((!b1.GetSourceAtom() && !b1.GetTargetAtom()));
  BOOST_TEST((b2.GetSourceAtom() == a1 && b2.GetTargetAtom() == a2));
  
}

BOOST_AUTO_TEST_CASE(aromaticity_get_set) {
  IXBond b = IXBond();
  BOOST_CHECK(!b.GetAromaticity());
  b.SetAromaticity(true);
  BOOST_CHECK(b.GetAromaticity());
  b.SetAromaticity(false);
  BOOST_CHECK(!b.GetAromaticity());
}

BOOST_AUTO_TEST_CASE(molecule_get_set) {
  IXBond b = IXBond();
  Molecule mol = Molecule(new IXMolecule());
  BOOST_CHECK(!b.GetMolecule());
  b.SetMolecule(mol);
  BOOST_CHECK(b.GetMolecule() == mol);
  mol.reset();
  BOOST_CHECK(!b.GetMolecule());
}

BOOST_AUTO_TEST_CASE(order_get_set) {
  using BondOrder = IXBond::Order;
  IXBond b = IXBond();
  BOOST_CHECK(b.GetOrder() == BondOrder::UNDEFINED);
  b.SetOrder(BondOrder::TRIPLE);
  BOOST_CHECK(b.GetOrder() == BondOrder::TRIPLE);
  b.SetOrder(BondOrder::TWOANDAHALF);
  BOOST_CHECK(b.GetOrder() == BondOrder::TWOANDAHALF);
}

BOOST_AUTO_TEST_CASE(atoms_get_set) {
  Atom a1 = Atom(new IXAtom());
  Atom a2 = Atom(new IXAtom());
  IXBond b = IXBond();
  BOOST_CHECK(b.NumAtoms() == 0);
  BOOST_CHECK(!b.GetSourceAtom());
  BOOST_CHECK(!b.GetTargetAtom());
  BOOST_CHECK(b.SetSourceAtom(a2));
  BOOST_CHECK(!b.SetTargetAtom(a2));
  BOOST_CHECK(b.SetTargetAtom(a1));
  BOOST_CHECK(!b.SetSourceAtom(a1));
  // Check correct assignment
  auto iters = b.GetAtomIters();
  BOOST_CHECK(std::distance(iters.first, iters.second) == 2);
  BOOST_CHECK(b.GetSourceAtom() == a2);
  BOOST_CHECK(b.GetTargetAtom() == a1);
  // Check correct swapping
  b.SwapSourceTarget();
  BOOST_CHECK(b.GetSourceAtom() == a1);
  BOOST_CHECK(b.GetTargetAtom() == a2);
  BOOST_CHECK(b.NumAtoms() == 2);
  a1.reset();
  BOOST_CHECK(b.NumAtoms() == 1);
  iters = b.GetAtomIters();
  BOOST_CHECK(std::distance(iters.first, iters.second) == 2);
}

BOOST_AUTO_TEST_CASE(stereochemistry_get_set) {
  using BondStereo = IXBond::Stereo;
  IXBond b = IXBond();
  BOOST_CHECK(b.GetStereochemistry() == BondStereo::UNDEFINED);
  b.SetStereochemistry(BondStereo::NONE);
  BOOST_CHECK(b.GetStereochemistry() == BondStereo::NONE);
  b.SetStereochemistry(BondStereo::Z);
  BOOST_CHECK(b.GetStereochemistry() == BondStereo::Z);
}

BOOST_AUTO_TEST_CASE(tag_get_set) {
  IXBond b = IXBond();
  BOOST_CHECK(b.GetTag() == 0);
  b.SetTag(12);
  BOOST_CHECK(b.GetTag() == 12);
  b.SetTag(123);
  BOOST_CHECK(b.GetTag() == 123);
}

BOOST_AUTO_TEST_CASE(angle_add_remove) {
  IXBond b = IXBond();
  Angle ang1 = Angle(new IXAngle());
  Angle ang2 = Angle(new IXAngle());
  Angle ang3 = Angle(new IXAngle());
  
  BOOST_CHECK(b.NumAngles() == 0);
  BOOST_CHECK(b.AddAngle(ang1));
  BOOST_CHECK(b.AddAngle(ang2));
  BOOST_CHECK(b.AddAngle(ang3));
  BOOST_CHECK(!b.AddAngle(Angle()));
  BOOST_CHECK(b.NumAngles() == 3);
  BOOST_CHECK(!b.AddAngle(ang3));
  BOOST_CHECK(b.NumAngles() == 3);
  BOOST_CHECK(b.RemoveAngle(ang3));
  BOOST_CHECK(!b.RemoveAngle(ang3));
  BOOST_CHECK(!b.RemoveAngle(Angle()));
  BOOST_CHECK(b.NumAngles() == 2);
  ang1.reset();
  BOOST_CHECK(b.NumAngles() == 1);
  auto iters = b.GetAngleIters();
  BOOST_CHECK(std::distance(iters.first, iters.second) == 2);
  b.Cleanup();
  iters = b.GetAngleIters();
  BOOST_CHECK(std::distance(iters.first, iters.second) == 1);
}

BOOST_AUTO_TEST_CASE(dihedral_add_remove) {
  IXBond b = IXBond();
  Dihedral ang1 = Dihedral(new IXDihedral());
  Dihedral ang2 = Dihedral(new IXDihedral());
  Dihedral ang3 = Dihedral(new IXDihedral());
  
  BOOST_CHECK(b.NumDihedrals() == 0);
  BOOST_CHECK(b.AddDihedral(ang1));
  BOOST_CHECK(b.AddDihedral(ang2));
  BOOST_CHECK(b.AddDihedral(ang3));
  BOOST_CHECK(!b.AddDihedral(Dihedral()));
  BOOST_CHECK(b.NumDihedrals() == 3);
  BOOST_CHECK(!b.AddDihedral(ang3));
  BOOST_CHECK(b.NumDihedrals() == 3);
  BOOST_CHECK(b.RemoveDihedral(ang3));
  BOOST_CHECK(!b.RemoveDihedral(ang3));
  BOOST_CHECK(!b.RemoveDihedral(Dihedral()));
  BOOST_CHECK(b.NumDihedrals() == 2);
  ang1.reset();
  BOOST_CHECK(b.NumDihedrals() == 1);
  auto iters = b.GetDihedralIters();
  BOOST_CHECK(std::distance(iters.first, iters.second) == 2);
  b.Cleanup();
  iters = b.GetDihedralIters();
  BOOST_CHECK(std::distance(iters.first, iters.second) == 1);
}

BOOST_AUTO_TEST_CASE(printing_methods) {
  Atom atm1 = Atom(new IXAtom());
  Atom atm2 = Atom(new IXAtom());
  Bond b = Bond(new IXBond(atm1, atm2));
  // Test ToString
  BOOST_CHECK(b->ToString() == "Bond(Atom(, XX), Atom(, XX))");
  atm1->SetName("Romeo"); atm1->SetElement("U");
  atm2->SetName("Foxtrot"); atm2->SetElement("F");
  b->SwapSourceTarget();
  BOOST_CHECK(b->ToString() == "Bond(Atom(Foxtrot, F), Atom(Romeo, U))");
  atm1.reset();
  BOOST_CHECK(b->ToString() == "Bond(MALFORMED)");
  // Test ostreaming
  std::stringstream expected, actual;
  expected << "Bond(Atom(" << atm2->GetUniqueID() << ", F), )";
  actual << b;
  BOOST_CHECK(expected.str() == actual.str());
  b.reset();
  expected.str("");
  actual.str("");
  actual << b;
  BOOST_CHECK(expected.str() == actual.str());
}

BOOST_AUTO_TEST_CASE(cleaning_methods) {
  Atom atm1 = Atom(new IXAtom());
  Atom atm2 = Atom(new IXAtom());
  Molecule mol = Molecule(new IXMolecule());
  IXBond b = IXBond(atm1, atm2);
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;
  angles.reserve(5); dihedrals.reserve(5);
  for (int i = 0; i < 5; ++i) {
    angles.emplace_back(new IXAngle());
    dihedrals.emplace_back(new IXDihedral());
    BOOST_CHECK(b.AddAngle(angles.back()));
    BOOST_CHECK(b.AddDihedral(dihedrals.back()));
  }
  b.SetMolecule(mol);
  b.SetAromaticity(true);
  b.SetOrder(indigox::IXBond::Order::QUADRUPLE);
  b.SetStereochemistry(indigox::IXBond::Stereo::NONE);
  b.SetTag(12);
  BOOST_CHECK(b.GetAromaticity());
  BOOST_CHECK(b.GetMolecule() == mol);
  BOOST_CHECK(b.GetOrder() == IXBond::Order::QUADRUPLE);
  BOOST_CHECK(b.GetSourceAtom() == atm1);
  BOOST_CHECK(b.GetStereochemistry() == IXBond::Stereo::NONE);
  BOOST_CHECK(b.GetTag() == 12);
  BOOST_CHECK(b.GetTargetAtom() == atm2);
  BOOST_CHECK(b.NumAngles() == 5);
  BOOST_CHECK(b.NumAtoms() == 2);
  BOOST_CHECK(b.NumDihedrals() == 5);
  b.Clear();
  BOOST_CHECK(!b.GetAromaticity());
  BOOST_CHECK(!b.GetMolecule());
  BOOST_CHECK(b.GetOrder() == IXBond::Order::UNDEFINED);
  BOOST_CHECK(!b.GetSourceAtom());
  BOOST_CHECK(b.GetStereochemistry() == IXBond::Stereo::UNDEFINED);
  BOOST_CHECK(!b.GetTag());
  BOOST_CHECK(!b.GetTargetAtom());
  BOOST_CHECK(b.NumAngles() == 0);
  BOOST_CHECK(b.NumAtoms() == 0);
  BOOST_CHECK(b.NumDihedrals() == 0);
}


BOOST_AUTO_TEST_SUITE_END();
