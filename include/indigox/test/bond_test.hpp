#ifndef INDIGOX_TEST_BOND_TEST_HPP
#define INDIGOX_TEST_BOND_TEST_HPP

#include "atom_test.hpp"

namespace indigox::test {
  struct TestBond {
    // Typedefs
    using IAtoms = IXBond::BondAtoms;
    using Atoms = std::pair<Atom, Atom>;
    
    indigox::Bond imp;
    
    // Private wrapping functions
    TestBond() = delete;
    TestBond(Atom a, Atom b, Molecule m) : imp(new IXBond(a,b,m)) { }
    void Clear() { imp->Clear(); }
    
    // Public wrapping members
    uid_ GetTag() const { return imp->GetTag(); }
    Molecule GetMolecule() const { return imp->GetMolecule(); }
    BondOrder GetOrder() const { return imp->GetOrder(); }
    Atom GetSourceAtom() const { return imp->GetSourceAtom(); }
    Atom GetTargetAtom() const { return imp->GetTargetAtom(); }
    bool GetAromaticity() const { return imp->GetAromaticity(); }
    BondStereo GetStereochemistry() const { return imp->GetStereochemistry(); }
    string_ ToString() const { return imp->ToString(); }
    void SetTag(uid_ t) { imp->SetTag(t); }
    void SetOrder(BondOrder o) { imp->SetOrder(o); }
    void SwapSourceTarget() { imp->SwapSourceTarget(); }
    void SetAromaticity(bool a) { imp->SetAromaticity(a); }
    void SetStereochemistry(BondStereo s) { imp->SetStereochemistry(s); }
    Atoms GetAtoms() const { return imp->GetAtoms(); }
    size_ NumAtoms() const { return imp->NumAtoms(); }
    uid_ GetUniqueID() const { return imp->GetUniqueID(); }
    
    // Internals access
    _Molecule get_mol() const { return imp->_mol; }
    uid_ get_tag() const { return imp->_tag; }
    BondOrder get_order() const { return imp->_order; }
    bool get_aromatic() const { return imp->_aromatic; }
    BondStereo get_stereo() const { return imp->_stereo; }
    const IAtoms& get_atms() const { return imp->_atms; }
  };
  
  inline TestBond CreateGenericTestBond() {
    return TestBond(CreateGenericTestAtom().imp,
                    CreateGenericTestAtom().imp,
                    Molecule());
  }
  
  struct BondTestFixture {
    Molecule mol = CreateMolecule();
    Atom a = CreateGenericTestAtom().imp;
    Atom b = CreateGenericTestAtom().imp;
    TestBond bnd = TestBond(a,b,mol);
    BondTestFixture() {
      a->SetTag(0); b->SetTag(1);
      a->SetElement("C"); b->SetElement("O");
    }
  };
}

#endif /* INDIGOX_TEST_BOND_TEST_HPP */
