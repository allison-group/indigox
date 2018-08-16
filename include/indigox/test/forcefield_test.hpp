#ifndef INDIGOX_TEST_FORCEFIELD_TEST_HPP
#define INDIGOX_TEST_FORCEFIELD_TEST_HPP

namespace indigox::test {
  
  struct TestFFDihedral {
    // Typedefs
    using DataStore = indigox::IXFFDihedral::DataStore;
    using AllowMask = indigox::IXFFDihedral::AllowedMask;
    using AllowEnum = indigox::IXFFDihedral::AllowEnum;
    using StoreEnum = indigox::IXFFDihedral::StoreEnum;
    
    indigox::FFDihedral imp;
    
    // Private wrapping functions
    TestFFDihedral() = delete;
    TestFFDihedral(DihedralType t, int_ i, std::initializer_list<float_> l)
    : imp(new IXFFDihedral(t,i,l)) { }
    
    // Public wrapping functions
    float_ GetPhaseShift() const { return imp->GetPhaseShift(); }
    float_ GetForceConstant() const { return imp->GetForceConstant(); }
    uint_ GetMultiplicity() const { return imp->GetMultiplicity(); }
    float_ GetIdealAngle() const { return imp->GetIdealAngle(); }
    DihedralType GetType() const { return imp->GetType(); }
    int_ GetID() const { return imp->GetID(); }
    
    // Internals access
    DihedralType get_type() const { return imp->_type; }
    int_ get_id() const { return imp->_id; }
    const DataStore& get_dat() const { return imp->_dat; }
    AllowMask& get_mask() const { return imp->_mask; }
    
  };
  
  inline TestFFDihedral CreateGenericTestFFDihedral() {
    return TestFFDihedral(DihedralType::Improper, 3, {0.95, 1.37});
  }
  
  struct TestFFAngle {
    // Typedefs
    using DataStore = indigox::IXFFAngle::DataStore;
    using AllowMask = indigox::IXFFAngle::AllowedMask;
    using AllowEnum = indigox::IXFFAngle::AllowEnum;
    using StoreEnum = indigox::IXFFAngle::StoreEnum;
    
    indigox::FFAngle imp;
    
    // Private wrapping functions
    TestFFAngle() = delete;
    TestFFAngle(AngleType t, int_ i, std::initializer_list<float_> l)
    : imp(new IXFFAngle(t,i,l)) { }
    
    // Public wrapping functions
    float_ GetForceConstant() const { return imp->GetForceConstant(); }
    float_ GetIdealAngle() const { return imp->GetIdealAngle(); }
    AngleType GetType() const { return imp->GetType(); }
    int_ GetID() const { return imp->GetID(); }
    
    // Internals access
    AngleType get_type() const { return imp->_type; }
    int_ get_id() const { return imp->_id; }
    const DataStore& get_dat() const { return imp->_dat; }
    AllowMask& get_mask() const { return imp->_mask; }
    
  };
  
  inline TestFFAngle CreateGenericTestFFAngle() {
    return TestFFAngle(AngleType::Harmonic, 3, {0.56,0.0098});
  }
  
  struct TestFFBond {
    // Typedefs
    using DataStore = indigox::IXFFBond::DataStore;
    using AllowMask = indigox::IXFFBond::AllowedMask;
    using AllowEnum = indigox::IXFFBond::AllowEnum;
    using StoreEnum = indigox::IXFFBond::StoreEnum;
    
    indigox::FFBond imp;
    
    // Private wrapping functions
    TestFFBond() = delete;
    TestFFBond(BondType t, int_ i, std::initializer_list<float_> l)
    : imp(new IXFFBond(t,i,l)) { }
    
    // Public wrapping functions
    float_ GetForceConstant() const { return imp->GetForceConstant(); }
    float_ GetIdealLength() const { return imp->GetIdealLength(); }
    BondType GetType() const { return imp->GetType(); }
    int_ GetID() const { return imp->GetID(); }
    
    // Internals access
    BondType get_type() const { return imp->_type; }
    int_ get_id() const { return imp->_id; }
    const DataStore& get_dat() const { return imp->_dat; }
    AllowMask& get_mask() const { return imp->_mask; }
  };
  
  inline TestFFBond CreateGenericTestFFBond() {
    return TestFFBond(BondType::Harmonic, 4, {0.567, 9.732});
  }
  
  struct TestFFAtom {
    // Typedefs
    
    indigox::FFAtom imp;
    
    // Private wrapping functions
    TestFFAtom() = delete;
    TestFFAtom(int_ i, string_ n) : imp(new IXFFAtom(i,n)) { }
    
    // Public wrapping functions
    int_ GetID() const { return imp->GetID(); }
    string_ GetName() const { return imp->GetName(); }
    
    // Internals access
    int_ get_id() const { return imp->_id; }
    string_ get_name() const { return imp->_name; }
  };
  
  inline TestFFAtom CreateGenericTestFFAtom() {
    return TestFFAtom(7, "GEN");
  }
  
  struct TestForcefield {
    // Typedefs
    using AtomTypes = indigox::IXForcefield::AtomTypes;
    using BondTypes = indigox::IXForcefield::BondTypes;
    using AngleTypes = indigox::IXForcefield::AngleTypes;
    using DihedralTypes = indigox::IXForcefield::DihedralTypes;
    using BndT = BondType;
    using AngT = AngleType;
    using DhdT = DihedralType;
    using ParmT = std::initializer_list<float_>;
    
    indigox::Forcefield imp;
    
    // Private wrapping functions
    FFBond NewBondType(BndT t, int_ i, ParmT p) {
      return imp->NewBondType(t, i, p);
    }
    FFAngle NewAngleType(AngT t, int_ i, ParmT p) {
      return imp->NewAngleType(t, i, p);
    }
    FFDihedral NewDihedralType(DhdT t, int_ i, ParmT p) {
      return imp->NewDihedralType(t, i, p);
    }
    
    // Public wrapping functions
    TestForcefield() = delete;
    TestForcefield(FFFamily f, string_ n) : imp(new IXForcefield(f,n)) { }
    void ReserveAtomTypes(size_ s) { imp->ReserveAtomTypes(s); }
    void ReserveBondTypes(BndT t, size_ s) { imp->ReserveBondTypes(t, s); }
    void ReserveAngleTypes(AngT t, size_ s) { imp->ReserveAngleTypes(t, s); }
    void ReserveDihedralTypes(DhdT t, size_ s) { imp->ReserveDihedralTypes(t, s); }
    FFFamily GetFamily() const { return imp->GetFamily(); }
    string_ GetName() const { return imp->GetName(); }
    FFAtom NewAtomType(int_ i, string_ n) { return imp->NewAtomType(i,n); }
    FFAtom GetAtomType(string_ n) const { return imp->GetAtomType(n); }
    FFAtom GetAtomType(int_ i) const { return imp->GetAtomType(i); }
    FFBond NewHarmonicBondType(int_ i, float_ k, float_ l) {
      return imp->NewHarmonicBondType(i, k, l);
    }
    FFBond NewQuarticBondType(int_ i, float_ k, float_ l) {
      return imp->NewQuarticBondType(i, k, l);
    }
    FFBond GetBondType(BndT t, int_ i) const { return imp->GetBondType(t,i); }
    FFBond GetHarmonicBondType(int_ i) const {
      return imp->GetHarmonicBondType(i);
    }
    FFBond GetQuarticBondType(int_ i) const {
      return imp->GetQuarticBondType(i);
    }
    FFAngle NewHarmonicAngleType(int_ i, float_ k, float_ t) {
      return imp->NewHarmonicAngleType(i, k, t);
    }
    FFAngle NewCosineHarmonicAngleType(int_ i, float_ k, float_ t) {
      return imp->NewCosineHarmonicAngleType(i, k, t);
    }
    FFAngle GetAngleType(AngT t, int_ i) const { return imp->GetAngleType(t, i); }
    FFAngle GetHarmonicAngleType(int_ i) const {
      return imp->GetHarmonicAngleType(i);
    }
    FFAngle GetCosineHarmonicAngleType(int_ i) const {
      return imp->GetCosineHarmonicAngleType(i);
    }
    FFDihedral NewProperDihedralType(int_ i, float_ k, float_ p, uint_ m) {
      return imp->NewProperDihedralType(i, k, p, m);
    }
    FFDihedral NewImproperDihedralType(int_ i, float_ k, float_ t) {
      return imp->NewImproperDihedralType(i, k, t);
    }
    FFDihedral GetDihedralType(DhdT t, int_ i) const {
      return imp->GetDihedralType(t, i);
    }
    FFDihedral GetImproperDihedralType(int_ i) const {
      return imp->GetImproperDihedralType(i);
    }
    FFDihedral GetProperDihedralType(int_ i) const {
      return imp->GetProperDihedralType(i);
    }
    size_ NumAtomTypes() const { return imp->NumAtomTypes(); }
    size_ NumBondTypes() const { return imp->NumBondTypes(); }
    size_ NumAngleTypes() const { return imp->NumAngleTypes(); }
    size_ NumDihedralTypes() const { return imp->NumDihedralTypes(); }
  };
  
}

#endif
