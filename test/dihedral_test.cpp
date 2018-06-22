#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

#include "class_test_wrappers.hpp"

namespace indigox::test {
  struct DihedralFixture {
    Molecule mol;
    Atom a1, a2, a3, a4, a5, a6, a7, a8;
    indigox::test::IXDihedral d1, d2;
    
    DihedralFixture() : mol(CreateMolecule()), a1(indigox::test::IXAtom::GetNewAtom()),
    a2(indigox::test::IXAtom::GetNewAtom()), a3(indigox::test::IXAtom::GetNewAtom()),
    a4(indigox::test::IXAtom::GetNewAtom()), a5(indigox::test::IXAtom::GetNewAtom()),
    a6(indigox::test::IXAtom::GetNewAtom()), a7(indigox::test::IXAtom::GetNewAtom()),
    a8(indigox::test::IXAtom::GetNewAtom()), d1(a1,a2,a3,a4,Molecule()), d2(a7,a8,a6,a5,mol)
    {}
  };
}

using namespace indigox;

BOOST_FIXTURE_TEST_SUITE(ixdihedral, indigox::test::DihedralFixture);

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK_NO_THROW(test::IXDihedral());
  BOOST_CHECK_NO_THROW(test::IXDihedral(a2,a6,a8,a7, mol));
  BOOST_CHECK(d1.GetUniqueId());
  BOOST_CHECK(d2.GetUniqueId());
  BOOST_CHECK(d2.GetUniqueId() != d1.GetUniqueId());
}

BOOST_AUTO_TEST_CASE(molecule_get) {
  BOOST_CHECK(d1.GetMolecule() == Molecule());
  BOOST_CHECK(d2.GetMolecule() == mol);
  mol.reset();
  BOOST_CHECK(!d2.GetMolecule());
  BOOST_CHECK(d2.GetMolecule() == Molecule());
}

BOOST_AUTO_TEST_CASE(tag_get_set) {
  BOOST_CHECK(d1.GetTag() == 0);
  BOOST_CHECK(d2.GetTag() == 0);
  d2.SetTag(34);
  BOOST_CHECK(d2.GetTag() == 34);
}

BOOST_AUTO_TEST_CASE(atoms_get_swap) {
  BOOST_CHECK(d1.NumAtoms() == 4);
  BOOST_CHECK(d2.NumAtoms() == 4);
  BOOST_CHECK(d1.GetAtoms().first == a1);
  BOOST_CHECK(d1.GetAtoms().second == a2);
  BOOST_CHECK(d1.GetAtoms().third == a3);
  BOOST_CHECK(d1.GetAtoms().fourth == a4);
  BOOST_CHECK(d2.GetAtoms().first == a7);
  BOOST_CHECK(d2.GetAtoms().second == a8);
  BOOST_CHECK(d2.GetAtoms().third == a6);
  BOOST_CHECK(d2.GetAtoms().fourth == a5);
  a1.reset();
  BOOST_CHECK(d1.GetAtoms().first == Atom());
  BOOST_CHECK(d1.NumAtoms() == 4);
  
  d2.SwapOrder();
  BOOST_CHECK(d2.GetAtoms().first == a5);
  BOOST_CHECK(d2.GetAtoms().second == a6);
  BOOST_CHECK(d2.GetAtoms().third == a8);
  BOOST_CHECK(d2.GetAtoms().fourth == a7);
}

BOOST_AUTO_TEST_CASE(printing_methods) {
  // Test ToString
  BOOST_CHECK(d1.ToString() == "Dihedral(Atom(, XX), Atom(, XX), Atom(, XX), Atom(, XX))");
  a3.reset();
  BOOST_CHECK(d1.ToString() == "Dihedral(MALFORMED)");
}

BOOST_AUTO_TEST_CASE(cleaning_methods) {
  d2.SetTag(72);
  d2.Clear();
  BOOST_CHECK(d2.GetAtoms().first == Atom());
  BOOST_CHECK(d2.GetAtoms().second == Atom());
  BOOST_CHECK(d2.GetAtoms().third == Atom());
  BOOST_CHECK(d2.GetAtoms().fourth == Atom());
  BOOST_CHECK(d2.GetMolecule() == Molecule());
  BOOST_CHECK(d2.GetTag() == 0);
  BOOST_CHECK(d2.NumAtoms() == 4);
}

BOOST_AUTO_TEST_SUITE_END();
