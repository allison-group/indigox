#ifndef INDIGOX_TEST_DIHEDRAL_TEST_HPP
#define INDIGOX_TEST_DIHEDRAL_TEST_HPP

#include "atom_test.hpp"

namespace indigox::test {
  struct TestDihedral {
//    // Typedefs
//    using IAtoms = IXDihedral::DihedAtoms;
//    using Atoms = stdx::quad<Atom, Atom, Atom, Atom>;
//
//    indigox::Dihedral imp;
//
//    // Private wrapping functions
//    TestDihedral() = delete;
//    TestDihedral(Atom i, Atom j, Atom k, Atom l, Molecule m)
//    : imp(new IXDihedral(i,j,k,l,m)) { }
//    void Clear() { imp->Clear(); }
//
//    // Public wrapping functions
//    uint32_t GetUniqueID() const { return imp->GetUniqueID(); }
//    uint32_t GetTag() const { return imp->GetTag(); }
//    Molecule GetMolecule() const { return imp->GetMolecule(); }
//    Atoms GetAtoms() const { return imp->GetAtoms(); }
//    size_t NumAtoms() const { return imp->NumAtoms(); }
//    void SwapOrder() { imp->SwapOrder(); }
//    std::string ToString() const { return imp->ToString(); }
//    void SetTag(uint32_t t) { imp->SetTag(t); }
//    size_t GetIndex() const { return imp->GetIndex(); }
//    FFDihedral GetType() const { return imp->GetType(); }
//    void SetType(FFDihedral t) { imp->SetType(t); }
//
//    // Internals access
//    _Molecule get_mol() const { return imp->_mol; }
//    uint32_t get_tag() const { return imp->_tag; }
//    const IAtoms& get_atms() const { return imp->_atms; }
//    FFDihedral get_type() const { return imp->_type; }
  };
  
//  inline TestDihedral CreateGenericTestDihedral() {
//    return TestDihedral(CreateGenericTestAtom().imp,
//                        CreateGenericTestAtom().imp,
//                        CreateGenericTestAtom().imp,
//                        CreateGenericTestAtom().imp,
//                        Molecule());
//  }
  
  struct DihedralTestFixture {
//    Molecule mol = CreateMolecule();
//    Atom a = CreateGenericTestAtom().imp;
//    Atom b = CreateGenericTestAtom().imp;
//    Atom c = CreateGenericTestAtom().imp;
//    Atom d = CreateGenericTestAtom().imp;
//    TestDihedral dhd = TestDihedral(a,b,c,d,mol);
//    FFDihedral fftype;
//    DihedralTestFixture();
  };
}

#endif /* INDIGOX_TEST_DIHEDRAL_TEST_HPP */
