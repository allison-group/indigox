#ifndef INDIGOX_TEST_ATOM_TEST_HPP
#define INDIGOX_TEST_ATOM_TEST_HPP

namespace indigox::test {
  struct TestAtom {
    indigox::IXAtom a;
    Atom imp;
    typedef indigox::IXAtom::AtomAngleIter AtmAngIter;
    typedef indigox::IXAtom::AtomBondIter AtmBndIter;
    typedef indigox::IXAtom::AtomDihedralIter AtmDhdIter;
    
    // private wrapping members
    TestAtom() = delete;
    TestAtom(Molecule m) : a(m), imp(new IXAtom(m)) { }
    inline void AddAngle(Angle i) { a.AddAngle(i); }
    inline void AddBond(Bond i) { a.AddBond(i); }
    inline void AddDihedral(Dihedral i) { a.AddDihedral(i); }
    inline void Clear() { a.Clear(); }
    inline void RemoveAngle(Angle i) { a.RemoveAngle(i); }
    inline void RemoveBond(Bond i) { a.RemoveBond(i); }
    inline void RemoveDihedral(Dihedral i) { a.RemoveDihedral(i); }
    
    // public wrapping members
    inline uint_ AddImplicitHydrogen() { return a.AddImplicitHydrogen(); }
    inline std::pair<AtmAngIter, AtmAngIter> GetAngleIters() { return a.GetAngleIters(); }
    inline bool GetAromaticity() { return a.GetAromaticity(); }
    inline std::pair<AtmBndIter, AtmBndIter> GetBondIters() { return a.GetBondIters(); }
    inline std::pair<AtmDhdIter, AtmDhdIter> GetDihedralIters() { return a.GetDihedralIters(); }
    inline Element GetElement() { return a.GetElement(); }
    inline int_ GetFormalCharge() { return a.GetFormalCharge(); }
    inline uint_ GetImplicitCount() { return a.GetImplicitCount(); }
    inline Molecule GetMolecule() { return a.GetMolecule(); }
    inline string_ GetName() { return a.GetName(); }
    inline float_ GetPartialCharge() { return a.GetPartialCharge(); }
    inline AtomStereo GetStereochemistry() { return a.GetStereochemistry(); }
    inline uint_ GetTag() { return a.GetTag(); }
    inline uid_ GetUniqueID() { return a.GetUniqueID(); }
    inline Vec3 GetVector() { return a.GetVector(); }
    inline float_ GetX() { return a.GetX(); }
    inline float_ GetY() { return a.GetY(); }
    inline float_ GetZ() { return a.GetZ(); }
    inline size_ NumAngles() { return a.NumAngles(); }
    inline size_ NumBonds() { return a.NumBonds(); }
    inline size_ NumDihedrals() { return a.NumDihedrals(); }
    inline uint_ RemoveImplicitHydrogen() { return a.RemoveImplicitHydrogen(); }
    inline void SetAromaticity(bool i) { a.SetAromaticity(i); }
    inline void SetElement(Element i) { a.SetElement(i); }
    inline void SetElement(string_ i) { a.SetElement(i); }
    inline void SetElement(uint_ i) { a.SetElement(i); }
    inline void SetFormalCharge(int_ i) { a.SetFormalCharge(i); }
    inline void SetImplicitCount(uint_ i) { a.SetImplicitCount(i); }
    inline void SetName(string_ i) { a.SetName(i); }
    inline void SetPartialCharge(float_ i) { a.SetPartialCharge(i); }
    inline void SetPosition(float_ x, float_ y, float_ z) { a.SetPosition(x,y,z); }
    inline void SetStereochemistry(AtomStereo i) { a.SetStereochemistry(i); }
    inline void SetTag(uint_ i) { a.SetTag(i); }
    inline void SetX(float_ i) { a.SetX(i); }
    inline void SetY(float_ i) { a.SetY(i); }
    inline void SetZ(float_ i) { a.SetZ(i); }
    inline string_ ToString() { return a.ToString(); }
    
  };
  
  inline TestAtom CreateGenericTestAtom() {
    return TestAtom(Molecule());
  }
}

#endif /* INDIGOX_TEST_ATOM_TEST_HPP */
