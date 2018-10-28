#ifndef INDIGOX_TEST_ANGLE_TEST_HPP
#define INDIGOX_TEST_ANGLE_TEST_HPP

#include "atom_test.hpp"

namespace indigox::test {
  struct TestAngle {
//    // Typedefs
//    using IAtoms = IXAngle::AngleAtoms;
//    using Atoms = stdx::triple<Atom, Atom, Atom>;
//
//    indigox::Angle imp;
//
//    // Private wrapping functions
//    TestAngle() = delete;
//    TestAngle(Atom i, Atom j, Atom k, Molecule m) : imp(new IXAngle(i,j,k,m)) { }
//    void Clear() { imp->Clear(); }
//
//    // Public wrapping members
//    Atoms GetAtoms() const { return imp->GetAtoms(); }
//    Molecule GetMolecule() const { return imp->GetMolecule(); }
//    uint32_t GetTag() const { return imp->GetTag(); }
//    uint32_t GetUniqueID() const { return imp->GetUniqueID(); }
//    size_t NumAtoms() const { return imp->NumAtoms(); }
//    void SetTag(uint32_t t) { imp->SetTag(t); }
//    void SwapOrder() { imp->SwapOrder(); }
//    std::string ToString() const { return imp->ToString(); }
//    size_t GetIndex() const { return imp->GetIndex(); }
//    FFAngle GetType() const { return imp->GetType(); }
//    void SetType(FFAngle t) { imp->SetType(t); }
//
//    // Internals access
//    _Molecule get_mol() const { return imp->_mol; }
//    uint32_t get_tag() const { return imp->_tag; }
//    const IAtoms& get_atms() const { return imp->_atms; }
//    FFAngle get_type() const { return imp->_type; }
  };
  
//  inline TestAngle CreateGenericTestAngle() {
//    return TestAngle(CreateGenericTestAtom().imp,
//                     CreateGenericTestAtom().imp,
//                     CreateGenericTestAtom().imp,
//                     Molecule());
//  }
  
  struct AngleTestFixture {
//    Molecule mol = CreateMolecule();
//    Atom a = CreateGenericTestAtom().imp;
//    Atom b = CreateGenericTestAtom().imp;
//    Atom c = CreateGenericTestAtom().imp;
//    TestAngle ang = TestAngle(a,b,c,mol);
//    AngleTestFixture() {
//      a->SetTag(0); b->SetTag(1); c->SetTag(2);
//      a->SetElement("C"); b->SetElement("O"); c->SetElement("N");
//    }
  };
}

#endif  // INDIGOX_TEST_ANGLE_TEST_HPP
