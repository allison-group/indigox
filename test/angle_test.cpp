#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

#include "class_test_wrappers.hpp"

namespace indigox::test {
  struct AngleFixture {
    Molecule mol;
    Atom a1, a2, a3, a4, a5, a6;
    indigox::test::IXAngle n1, n2;
    
    AngleFixture() : mol(CreateMolecule()), a1(indigox::test::IXAtom::GetNewAtom()),
    a2(indigox::test::IXAtom::GetNewAtom()), a3(indigox::test::IXAtom::GetNewAtom()),
    a4(indigox::test::IXAtom::GetNewAtom()), a5(indigox::test::IXAtom::GetNewAtom()),
    a6(indigox::test::IXAtom::GetNewAtom()), n1(a1,a2,a3,Molecule()), n2(a6,a5,a4,mol) {
    }
  };
}

using namespace indigox;

BOOST_FIXTURE_TEST_SUITE(ixangle, indigox::test::AngleFixture);

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK_NO_THROW(test::IXAngle());
  BOOST_CHECK_NO_THROW(test::IXAngle(a2,a3,a4,mol));
  BOOST_CHECK(n1.GetUniqueID());
  BOOST_CHECK(n2.GetUniqueID());
  BOOST_CHECK(n1.GetUniqueID() != n2.GetUniqueID());
}

BOOST_AUTO_TEST_CASE(molecule_get) {
  BOOST_CHECK(n1.GetMolecule() == Molecule());
  BOOST_CHECK(n2.GetMolecule() == mol);
  mol.reset();
  BOOST_CHECK(!n2.GetMolecule());
  BOOST_CHECK(n2.GetMolecule() == Molecule());
}

BOOST_AUTO_TEST_CASE(tag_get_set) {
  BOOST_CHECK(n1.GetTag() == 0);
  BOOST_CHECK(n2.GetTag() == 0);
  n2.SetTag(12);
  BOOST_CHECK(n2.GetTag() == 12);
}

BOOST_AUTO_TEST_CASE(atoms_get_swap) {
  BOOST_CHECK(n1.NumAtoms() == 3);
  BOOST_CHECK(n2.NumAtoms() == 3);
  BOOST_CHECK(n1.GetAtoms().first == a1);
  BOOST_CHECK(n2.GetAtoms().first == a6);
  BOOST_CHECK(n1.GetAtoms().second == a2);
  BOOST_CHECK(n2.GetAtoms().second == a5);
  BOOST_CHECK(n1.GetAtoms().third == a3);
  BOOST_CHECK(n2.GetAtoms().third == a4);
  a1.reset();
  BOOST_CHECK(n1.GetAtoms().first == Atom());
  BOOST_CHECK(n1.NumAtoms() == 3);
  
  n2.SwapOrder();
  BOOST_CHECK(n2.GetAtoms().first == a4);
  BOOST_CHECK(n2.GetAtoms().second == a5);
  BOOST_CHECK(n2.GetAtoms().third == a6);
}

BOOST_AUTO_TEST_CASE(printing_methods) {
  // Test ToString
  BOOST_CHECK(n1.ToString() == "Angle(Atom(, XX), Atom(, XX), Atom(, XX))");
  a1->SetName("A1"); a1->SetElement("B");
  a2->SetName("A2"); a2->SetElement("W");
  a3->SetName("A3");
  BOOST_CHECK(n1.ToString() == "Angle(Atom(A1, B), Atom(A2, W), Atom(A3, XX))");
  a3.reset();
  BOOST_CHECK(n1.ToString() == "Angle(MALFORMED)");
}

BOOST_AUTO_TEST_CASE(cleaning_methods) {
  n2.SetTag(72);
  n2.Clear();
  BOOST_CHECK(n2.GetAtoms().first == Atom());
  BOOST_CHECK(n2.GetAtoms().second == Atom());
  BOOST_CHECK(n2.GetAtoms().third == Atom());
  BOOST_CHECK(n2.GetMolecule() == Molecule());
  BOOST_CHECK(n2.GetTag() == 0);
  BOOST_CHECK(n2.NumAtoms() == 3);
}

BOOST_AUTO_TEST_SUITE_END();
