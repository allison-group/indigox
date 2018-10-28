#ifndef INDIGOX_TEST_BOND_TEST_HPP
#define INDIGOX_TEST_BOND_TEST_HPP

#include "atom_test.hpp"

namespace indigox::test {
  struct TestBond {
//    // Typedefs
//    using IAtoms = IXBond::BondAtoms;
//    using Atoms = std::pair<Atom, Atom>;
//
//    indigox::Bond imp;
//
//    // Private wrapping functions
//    TestBond() = delete;
//    TestBond(Atom a, Atom b, Molecule m) : imp(new IXBond(a,b,m)) { }
//    void Clear() { imp->Clear(); }
//
//    // Public wrapping members
//    uint32_t GetTag() const { return imp->GetTag(); }
//    Molecule GetMolecule() const { return imp->GetMolecule(); }
//    BondOrder GetOrder() const { return imp->GetOrder(); }
//    Atom GetSourceAtom() const { return imp->GetSourceAtom(); }
//    Atom GetTargetAtom() const { return imp->GetTargetAtom(); }
//    bool GetAromaticity() const { return imp->GetAromaticity(); }
//    BondStereo GetStereochemistry() const { return imp->GetStereochemistry(); }
//    std::string ToString() const { return imp->ToString(); }
//    void SetTag(uint32_t t) { imp->SetTag(t); }
//    void SetOrder(BondOrder o) { imp->SetOrder(o); }
//    void SwapSourceTarget() { imp->SwapSourceTarget(); }
//    void SetAromaticity(bool a) { imp->SetAromaticity(a); }
//    void SetStereochemistry(BondStereo s) { imp->SetStereochemistry(s); }
//    Atoms GetAtoms() const { return imp->GetAtoms(); }
//    size_t NumAtoms() const { return imp->NumAtoms(); }
//    uint32_t GetUniqueID() const { return imp->GetUniqueID(); }
//    size_t GetIndex() const { return imp->GetIndex(); }
//    FFBond GetType() const { return imp->GetType(); }
//    void SetType(FFBond t) { imp->SetType(t); }
//
//    // Internals access
//    _Molecule get_mol() const { return imp->_mol; }
//    uint32_t get_tag() const { return imp->_tag; }
//    BondOrder get_order() const { return imp->_order; }
//    bool get_aromatic() const { return imp->_aromatic; }
//    BondStereo get_stereo() const { return imp->_stereo; }
//    const IAtoms& get_atms() const { return imp->_atms; }
//    FFBond get_type() const { return imp->_type; }
  };
  
//  inline TestBond CreateGenericTestBond() {
//    return TestBond(CreateGenericTestAtom().imp,
//                    CreateGenericTestAtom().imp,
//                    Molecule());
//  }
//
//  inline TestBond CreateGenericTestBond(const Atom& u, const Atom& v) {
//    return TestBond(u,v,Molecule());
//  }
  
  struct BondTestFixture {
//    Molecule mol = CreateMolecule();
//    Atom a = CreateGenericTestAtom().imp;
//    Atom b = CreateGenericTestAtom().imp;
//    TestBond bnd = TestBond(a,b,mol);
//    FFBond fftype;
//    BondTestFixture();
  };
}

#endif /* INDIGOX_TEST_BOND_TEST_HPP */
