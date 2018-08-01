#ifndef INDIGOX_TEST_DIHEDRAL_TEST_HPP
#define INDIGOX_TEST_DIHEDRAL_TEST_HPP

#include "atom_test.hpp"

namespace indigox::test {
  struct TestDihedral {
    // Typedefs
    using IAtoms = IXDihedral::DihedAtoms;
    using Atoms = stdx::quad<Atom, Atom, Atom, Atom>;
    
    indigox::Dihedral imp;
    
    // Private wrapping functions
    TestDihedral() = delete;
    TestDihedral(Atom i, Atom j, Atom k, Atom l, Molecule m)
    : imp(new IXDihedral(i,j,k,l,m)) { }
    void Clear() { imp->Clear(); }
    
    // Public wrapping functions
    uid_ GetUniqueID() const { return imp->GetUniqueID(); }
    uid_ GetTag() const { return imp->GetTag(); }
    Molecule GetMolecule() const { return imp->GetMolecule(); }
    Atoms GetAtoms() const { return imp->GetAtoms(); }
    size_ NumAtoms() const { return imp->NumAtoms(); }
    void SwapOrder() { imp->SwapOrder(); }
    string_ ToString() const { return imp->ToString(); }
    void SetTag(uid_ t) { imp->SetTag(t); }
    size_ GetIndex() const { return imp->GetIndex(); }
    
    // Internals access
    _Molecule get_mol() const { return imp->_mol; }
    uid_ get_tag() const { return imp->_tag; }
    const IAtoms& get_atms() const { return imp->_atms; }
  };
  
  inline TestDihedral CreateGenericTestDihedral() {
    return TestDihedral(CreateGenericTestAtom().imp,
                        CreateGenericTestAtom().imp,
                        CreateGenericTestAtom().imp,
                        CreateGenericTestAtom().imp,
                        Molecule());
  }
  
  struct DihedralTestFixture {
    Molecule mol = CreateMolecule();
    Atom a = CreateGenericTestAtom().imp;
    Atom b = CreateGenericTestAtom().imp;
    Atom c = CreateGenericTestAtom().imp;
    Atom d = CreateGenericTestAtom().imp;
    TestDihedral dhd = TestDihedral(a,b,c,d,mol);
    DihedralTestFixture() {
      a->SetTag(0); b->SetTag(1); c->SetTag(2); d->SetTag(3);
      a->SetElement("C"); b->SetElement("O");
      c->SetElement("F"); d->SetElement("N");
    }
    
  };
}

#endif /* INDIGOX_TEST_DIHEDRAL_TEST_HPP */
