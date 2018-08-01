#ifndef INDIGOX_TEST_ATOM_TEST_HPP
#define INDIGOX_TEST_ATOM_TEST_HPP

#include <random>

namespace indigox::test {
  struct TestAtom {
    // Typedefs
    using IAngles = IXAtom::AtomAngles;
    using IBonds = IXAtom::AtomBonds;
    using IDihedrals = IXAtom::AtomDihedrals;
    using Angles = std::pair<IXAtom::AtomAngleIter, IXAtom::AtomAngleIter>;
    using Bonds = std::pair<IXAtom::AtomBondIter, IXAtom::AtomBondIter>;
    using Dihedrals = std::pair<IXAtom::AtomDihedralIter, IXAtom::AtomDihedralIter>;
    
    indigox::Atom imp;
    
    // Private wrapping members
    TestAtom() = delete;
    TestAtom(Molecule m) : imp(new IXAtom(m)) { }
    void AddAngle(Angle a) { imp->AddAngle(a); }
    void AddBond(Bond b) { imp->AddBond(b); }
    void AddDihedral(Dihedral d) { imp->AddDihedral(d); }
    void Clear() { imp->Clear(); }
    void RemoveAngle(Angle a) { imp->RemoveAngle(a); }
    void RemoveBond(Bond b) { imp->RemoveBond(b); }
    void RemoveDihedral(Dihedral d) { imp->RemoveDihedral(d); }
    
    // Public wrapping members
    uid_ GetUniqueID() const { return imp->GetUniqueID(); }
    Element GetElement() const { return imp->GetElement(); }
    int_ GetFormalCharge() const { return imp->GetFormalCharge(); }
    float_ GetPartialCharge() const { return imp->GetPartialCharge(); }
    uint_ GetTag() const { return imp->GetTag(); }
    uint_ GetImplicitCount() const { return imp->GetImplicitCount(); }
    uint_ AddImplicitHydrogen() { return imp->AddImplicitHydrogen(); }
    uint_ RemoveImplicitHydrogen() { return imp->RemoveImplicitHydrogen(); }
    Molecule GetMolecule() const { return imp->GetMolecule(); }
    string_ GetName() const { return imp->GetName(); }
    float_ GetX() const { return imp->GetX(); }
    float_ GetY() const { return imp->GetY(); }
    float_ GetZ() const { return imp->GetZ(); }
    const Vec3& GetVector() const { return imp->GetVector(); }
    string_ ToString() const { return imp->ToString(); }
    void SetElement(Element e) { imp->SetElement(e); }
    void SetElement(string_ e) { imp->SetElement(e); }
    void SetElement(uint_ e) { imp->SetElement(e); }
    void SetFormalCharge(int_ q) { imp->SetFormalCharge(q); }
    void SetPartialCharge(float_ q) { imp->SetPartialCharge(q); }
    void SetImplicitCount(uint_ h) { imp->SetImplicitCount(h); }
    void SetTag(uint_ t) { imp->SetTag(t); }
    void SetName(string_ n) { imp->SetName(n); }
    void SetX(float_ x) { imp->SetX(x); }
    void SetY(float_ y) { imp->SetY(y); }
    void SetZ(float_ z) { imp->SetZ(z); }
    void SetPosition(float_ x, float_ y, float_ z) { imp->SetPosition(x,y,z); }
    void SetStereochemistry(AtomStereo s) { imp->SetStereochemistry(s); }
    void SetAromaticity(bool a) { imp->SetAromaticity(a); }
    AtomStereo GetStereochemistry() const { return imp->GetStereochemistry(); }
    bool GetAromaticity() const { return imp->GetAromaticity(); }
    Bonds GetBondIters() const { return imp->GetBondIters(); }
    Angles GetAngleIters() const { return imp->GetAngleIters(); }
    Dihedrals GetDihedralIters() const { return imp->GetDihedralIters(); }
    size_ NumBonds() const { return imp->NumBonds(); }
    size_ NumAngles() const { return imp->NumAngles(); }
    size_ NumDihedrals() const { return imp->NumDihedrals(); }
    size_ GetIndex() const { return imp->GetIndex(); }

    // Internals access
    _Molecule get_mol() const { return imp->_mol; }
    _Element get_elem() const { return imp->_elem; }
    int_ get_fc() const { return imp->_fc; }
    uint_ get_tag() const { return imp->_tag; }
    uint_ get_implicitH() const { return imp->_implicitH; }
    string_ get_name() const { return imp->_name; }
    const Vec3& get_pos() const { return imp->_pos; }
    float_ get_partial() const { return imp->_partial; }
    AtomStereo get_stereo() const { return imp->_stereo; }
    bool get_aromatic() const { return imp->_aromatic; }
    const IBonds& get_bnds() const { return imp->_bnds; }
    const IAngles& get_angs() const { return imp->_angs; }
    const IDihedrals& get_dhds() const { return imp->_dhds; }
  };
  
  inline TestAtom CreateGenericTestAtom() {
    return TestAtom(Molecule());
  }
  
  struct AtomTestFixture {
    Molecule mol = CreateMolecule();
    TestAtom atm;
    AtomTestFixture() : mol(CreateMolecule()), atm(mol) { }
  };
}

#endif /* INDIGOX_TEST_ATOM_TEST_HPP */
