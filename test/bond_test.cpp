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
    
    BondFixture() : mol(CreateMolecule()), a1(indigox::test::IXAtom::GetNewAtom()),
    a2(indigox::test::IXAtom::GetNewAtom(mol)), a3(indigox::test::IXAtom::GetNewAtom()),
    a4(indigox::test::IXAtom::GetNewAtom(mol)), b1(a1,a3,Molecule()), b2(a2,a4,mol) {
    }
    
  };
}

using namespace indigox;

BOOST_FIXTURE_TEST_SUITE(ixbond, indigox::test::BondFixture);

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK_NO_THROW(test::IXBond());
  BOOST_CHECK_NO_THROW(test::IXBond(a2,a3,mol));
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

BOOST_AUTO_TEST_CASE(printing_methods) {
  // Test ostream operator
  boost::test_tools::output_test_stream actual_output;
  actual_output << indigox::test::IXBond::GetNewBond();
  BOOST_CHECK(actual_output.is_equal("Bond(, )"));
  
  // Test ToString
  BOOST_CHECK(b1.ToString() == "Bond(Atom(, XX), Atom(, XX))");
  a3.reset();  // malformed checking
  BOOST_CHECK(b1.ToString() == "Bond(MALFORMED)");
}

BOOST_AUTO_TEST_CASE(cleaning_methods) {
  b2.SetOrder(BondOrder::TRIPLE);
  b2.SetStereochemistry(BondStereo::Z);
  b2.SetTag(72);
  b2.Clear();
  BOOST_CHECK(b2.GetAromaticity() == false);
  BOOST_CHECK(b2.GetMolecule() == Molecule());
  BOOST_CHECK(b2.GetOrder() == BondOrder::UNDEFINED);
  BOOST_CHECK(b2.GetSourceAtom() == Atom());
  BOOST_CHECK(b2.GetStereochemistry() == BondStereo::UNDEFINED);
  BOOST_CHECK(b2.GetTag() == 0);
  BOOST_CHECK(b2.GetTargetAtom() == Atom());
  BOOST_CHECK(b2.NumAtoms() == 2);
}


BOOST_AUTO_TEST_SUITE_END();
