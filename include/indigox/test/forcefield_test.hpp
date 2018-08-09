#ifndef INDIGOX_TEST_FORCEFIELD_TEST_HPP
#define INDIGOX_TEST_FORCEFIELD_TEST_HPP

namespace indigox::test {
  
  struct TestDihedralType {
    // Typedefs
    using DataStore = indigox::IXFFDihedralType::DataStore;
    using AllowMask = indigox::IXFFDihedralType::AllowedMask;
    using AllowEnum = indigox::IXFFDihedralType::AllowEnum;
    using StoreEnum = indigox::IXFFDihedralType::StoreEnum;
    
    indigox::FFDihedralType imp;
    
    // Private wrapping functions
    TestDihedralType() = delete;
    TestDihedralType(DihedralFunctionType t, int_ i, std::initializer_list<float_> l)
    : imp(new IXFFDihedralType(t,i,l)) { }
    
    // Public wrapping functions
    float_ GetPhaseShift() const { return imp->GetPhaseShift(); }
    float_ GetForceConstant() const { return imp->GetForceConstant(); }
    uint_ GetMultiplicity() const { return imp->GetMultiplicity(); }
    float_ GetIdealAngle() const { return imp->GetIdealAngle(); }
    DihedralFunctionType GetType() const { return imp->GetType(); }
    int_ GetID() const { return imp->GetID(); }
    
    // Internals access
    DihedralFunctionType get_type() const { return imp->_type; }
    int_ get_id() const { return imp->_id; }
    const DataStore& get_dat() const { return imp->_dat; }
    AllowMask& get_mask() const { return imp->_mask; }
    
  };
  
  struct TestAngleType {
    // Typedefs
    using DataStore = indigox::IXFFAngleType::DataStore;
    using AllowMask = indigox::IXFFAngleType::AllowedMask;
    using AllowEnum = indigox::IXFFAngleType::AllowEnum;
    using StoreEnum = indigox::IXFFAngleType::StoreEnum;
    
    indigox::FFAngleType imp;
    
    // Private wrapping functions
    TestAngleType() = delete;
    TestAngleType(AngleFunctionType t, int_ i, std::initializer_list<float_> l)
    : imp(new IXFFAngleType(t,i,l)) { }
    
    // Public wrapping functions
    float_ GetForceConstant() const { return imp->GetForceConstant(); }
    float_ GetIdealAngle() const { return imp->GetIdealAngle(); }
    AngleFunctionType GetType() const { return imp->GetType(); }
    int_ GetID() const { return imp->GetID(); }
    
    // Internals access
    AngleFunctionType get_type() const { return imp->_type; }
    int_ get_id() const { return imp->_id; }
    const DataStore& get_dat() const { return imp->_dat; }
    AllowMask& get_mask() const { return imp->_mask; }
    
  };
  
  struct TestBondType {
    // Typedefs
    using DataStore = indigox::IXFFBondType::DataStore;
    using AllowMask = indigox::IXFFBondType::AllowedMask;
    using AllowEnum = indigox::IXFFBondType::AllowEnum;
    using StoreEnum = indigox::IXFFBondType::StoreEnum;
    
    indigox::FFBondType imp;
    
    // Private wrapping functions
    TestBondType() = delete;
    TestBondType(BondFunctionType t, int_ i, std::initializer_list<float_> l)
    : imp(new IXFFBondType(t,i,l)) { }
    
    // Public wrapping functions
    float_ GetForceConstant() const { return imp->GetForceConstant(); }
    float_ GetIdealLength() const { return imp->GetIdealLength(); }
    BondFunctionType GetType() const { return imp->GetType(); }
    int_ GetID() const { return imp->GetID(); }
    
    // Internals access
    BondFunctionType get_type() const { return imp->_type; }
    int_ get_id() const { return imp->_id; }
    const DataStore& get_dat() const { return imp->_dat; }
    AllowMask& get_mask() const { return imp->_mask; }
  };
  
  struct TestAtomType {
    // Typedefs
    
    indigox::FFAtomType imp;
    
    // Private wrapping functions
    TestAtomType() = delete;
    TestAtomType(int_ i, string_ n) : imp(new IXFFAtomType(i,n)) { }
    
    // Public wrapping functions
    int_ GetID() const { return imp->GetID(); }
    string_ GetName() const { return imp->GetName(); }
    
    // Internals access
    int_ get_id() const { return imp->_id; }
    string_ get_name() const { return imp->_name; }
  };
  
  struct TestForcefield {
    // Typedefs
    using AtomTypes = indigox::IXForcefield::AtomTypes;
    using BondTypes = indigox::IXForcefield::BondTypes;
    using AngleTypes = indigox::IXForcefield::AngleTypes;
    using DihedralTypes = indigox::IXForcefield::DihedralTypes;
    using BndT = BondFunctionType;
    using AngT = AngleFunctionType;
    using DhdT = DihedralFunctionType;
    using ParmT = std::initializer_list<float_>;
    
    indigox::Forcefield imp;
    
    // Private wrapping functions
    FFBondType NewBondType(BndT t, int_ i, ParmT p) {
      return imp->NewBondType(t, i, p);
    }
    FFAngleType NewAngleType(AngT t, int_ i, ParmT p) {
      return imp->NewAngleType(t, i, p);
    }
    FFDihedralType NewDihedralType(DhdT t, int_ i, ParmT p) {
      return imp->NewDihedralType(t, i, p);
    }
    
    // Public wrapping functions
    TestForcefield() = delete;
    TestForcefield(ForcefieldFamily f, string_ n) : imp(new IXForcefield(f,n)) { }
    void ReserveAtomTypes(size_ s) { imp->ReserveAtomTypes(s); }
    void ReserveBondTypes(BndT t, size_ s) { imp->ReserveBondTypes(t, s); }
    void ReserveAngleTypes(AngT t, size_ s) { imp->ReserveAngleTypes(t, s); }
    void ReserveDihedralTypes(DhdT t, size_ s) { imp->ReserveDihedralTypes(t, s); }
    ForcefieldFamily GetFamily() const { return imp->GetFamily(); }
    string_ GetName() const { return imp->GetName(); }
    FFAtomType NewAtomType(int_ i, string_ n) { return imp->NewAtomType(i,n); }
    FFAtomType GetAtomType(string_ n) const { return imp->GetAtomType(n); }
    FFAtomType GetAtomType(int_ i) const { return imp->GetAtomType(i); }
    FFBondType NewHarmonicBondType(int_ i, float_ k, float_ l) {
      return imp->NewHarmonicBondType(i, k, l);
    }
    FFBondType NewQuarticBondType(int_ i, float_ k, float_ l) {
      return imp->NewQuarticBondType(i, k, l);
    }
    FFBondType GetBondType(BndT t, int_ i) const { return imp->GetBondType(t,i); }
    FFBondType GetHarmonicBondType(int_ i) const {
      return imp->GetHarmonicBondType(i);
    }
    FFBondType GetQuarticBondType(int_ i) const {
      return imp->GetQuarticBondType(i);
    }
    FFAngleType NewHarmonicAngleType(int_ i, float_ k, float_ t) {
      return imp->NewHarmonicAngleType(i, k, t);
    }
    FFAngleType NewCosineHarmonicAngleType(int_ i, float_ k, float_ t) {
      return imp->NewCosineHarmonicAngleType(i, k, t);
    }
    FFAngleType GetAngleType(AngT t, int_ i) const { return imp->GetAngleType(t, i); }
    FFAngleType GetHarmonicAngleType(int_ i) const {
      return imp->GetHarmonicAngleType(i);
    }
    FFAngleType GetCosineHarmonicAngleType(int_ i) const {
      return imp->GetCosineHarmonicAngleType(i);
    }
    FFDihedralType NewProperDihedralType(int_ i, float_ k, float_ p, uint_ m) {
      return imp->NewProperDihedralType(i, k, p, m);
    }
    FFDihedralType NewImproperDihedralType(int_ i, float_ k, float_ t) {
      return imp->NewImproperDihedralType(i, k, t);
    }
    FFDihedralType GetDihedralType(DhdT t, int_ i) const {
      return imp->GetDihedralType(t, i);
    }
    FFDihedralType GetImproperDihedralType(int_ i) const {
      return imp->GetImproperDihedralType(i);
    }
    FFDihedralType GetProperDihedralType(int_ i) const {
      return imp->GetProperDihedralType(i);
    }
    size_ NumAtomTypes() const { return imp->NumAtomTypes(); }
    size_ NumBondTypes() const { return imp->NumBondTypes(); }
    size_ NumAngleTypes() const { return imp->NumAngleTypes(); }
    size_ NumDihedralTypes() const { return imp->NumDihedralTypes(); }
  };
  
}

#endif
