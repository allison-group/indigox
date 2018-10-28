#ifndef INDIGOX_TEST_FORCEFIELD_TEST_HPP
#define INDIGOX_TEST_FORCEFIELD_TEST_HPP

namespace indigox::test {
  
  struct TestFFDihedral {
    // Typedefs
//    using DataStore = indigox::IXFFDihedral::DataStore;
//    using AllowMask = indigox::IXFFDihedral::AllowedMask;
//    using AllowEnum = indigox::IXFFDihedral::AllowEnum;
//    using StoreEnum = indigox::IXFFDihedral::StoreEnum;
//
//    indigox::FFDihedral imp;
//
//    // Private wrapping functions
//    TestFFDihedral() = delete;
//    TestFFDihedral(DihedralType t, int i, FFParam l, const Forcefield& f)
//    : imp(new IXFFDihedral(t,i,l,f)) { }
//    TestFFDihedral(DihedralType t, const Forcefield& f)
//    : imp(new IXFFDihedral(t, f)) { }
//
//    // Public wrapping functions
//    double GetPhaseShift() const { return imp->GetPhaseShift(); }
//    double GetForceConstant() const { return imp->GetForceConstant(); }
//    uint32_t GetMultiplicity() const { return imp->GetMultiplicity(); }
//    double GetIdealAngle() const { return imp->GetIdealAngle(); }
//    DihedralType GetType() const { return imp->GetType(); }
//    int GetID() const { return imp->GetID(); }
//    Forcefield GetForcefield() const { return imp->GetForcefield(); }
//
//    // Internals access
//    DihedralType get_type() const { return imp->_type; }
//    int get_id() const { return imp->_id; }
//    const DataStore& get_dat() const { return imp->_dat; }
//    AllowMask& get_mask() const { return imp->_mask; }
//    _Forcefield get_ff() const { return imp->_ff; }
    
  };
  
  struct TestFFAngle {
    // Typedefs
//    using DataStore = indigox::IXFFAngle::DataStore;
//    using AllowMask = indigox::IXFFAngle::AllowedMask;
//    using AllowEnum = indigox::IXFFAngle::AllowEnum;
//    using StoreEnum = indigox::IXFFAngle::StoreEnum;
//
//    indigox::FFAngle imp;
//
//    // Private wrapping functions
//    TestFFAngle() = delete;
//    TestFFAngle(AngleType t, int i, FFParam l, const Forcefield& f)
//    : imp(new IXFFAngle(t,i,l,f)) { }
//    TestFFAngle(AngleType t, const Forcefield& f) : imp(new IXFFAngle(t,f)) { }
//
//    // Public wrapping functions
//    double GetForceConstant() const { return imp->GetForceConstant(); }
//    double GetIdealAngle() const { return imp->GetIdealAngle(); }
//    AngleType GetType() const { return imp->GetType(); }
//    int GetID() const { return imp->GetID(); }
//    FFAngle GetLinkedType() const { return imp->GetLinkedType(); }
//    Forcefield GetForcefield() const { return imp->GetForcefield(); }
//
//    // Internals access
//    AngleType get_type() const { return imp->_type; }
//    int get_id() const { return imp->_id; }
//    const DataStore& get_dat() const { return imp->_dat; }
//    AllowMask& get_mask() const { return imp->_mask; }
//    FFAngle get_link() const { return imp->_link; }
//    void set_link(FFAngle a) { imp->_link = a; }
//    _Forcefield get_ff() const { return imp->_ff; }
    
  };
  
  struct TestFFBond {
    // Typedefs
//    using DataStore = indigox::IXFFBond::DataStore;
//    using AllowMask = indigox::IXFFBond::AllowedMask;
//    using AllowEnum = indigox::IXFFBond::AllowEnum;
//    using StoreEnum = indigox::IXFFBond::StoreEnum;
//
//    indigox::FFBond imp;
//
//    // Private wrapping functions
//    TestFFBond() = delete;
//    TestFFBond(BondType t, int i, FFParam l, const Forcefield& f)
//    : imp(new IXFFBond(t,i,l,f)) { }
//    TestFFBond(BondType t, const Forcefield& f) : imp(new IXFFBond(t,f)) { }
//
//    // Public wrapping functions
//    double GetForceConstant() const { return imp->GetForceConstant(); }
//    double GetIdealLength() const { return imp->GetIdealLength(); }
//    BondType GetType() const { return imp->GetType(); }
//    int GetID() const { return imp->GetID(); }
//    FFBond GetLinkedType() const { return imp->GetLinkedType(); }
//    Forcefield GetForcefield() const { return imp->GetForcefield(); }
//
//    // Internals access
//    BondType get_type() const { return imp->_type; }
//    int get_id() const { return imp->_id; }
//    const DataStore& get_dat() const { return imp->_dat; }
//    AllowMask& get_mask() const { return imp->_mask; }
//    FFBond get_link() const { return imp->_link; }
//    void set_link(FFBond a) { imp->_link = a; }
//    _Forcefield get_ff() const { return imp->_ff; }
  };
  
  struct TestFFAtom {
    // Typedefs
    
//    indigox::FFAtom imp;
//
//    // Private wrapping functions
//    TestFFAtom() = delete;
//    TestFFAtom(int i, std::string n, const Forcefield& f)
//    : imp(new IXFFAtom(i,n,f)) { }
//
//    // Public wrapping functions
//    int GetID() const { return imp->GetID(); }
//    std::string GetName() const { return imp->GetName(); }
//    Forcefield GetForcefield() const { return imp->GetForcefield(); }
//
//    // Internals access
//    int get_id() const { return imp->_id; }
//    std::string get_name() const { return imp->_name; }
//    _Forcefield get_ff() const { return imp->_ff; }
  };
  
  struct TestForcefield {
    // Typedefs
//    using AtmTypes = indigox::IXForcefield::AtomTypes;
//    using BndTypes = indigox::IXForcefield::BondTypes;
//    using AngTypes = indigox::IXForcefield::AngleTypes;
//    using DhdTypes = indigox::IXForcefield::DihedralTypes;
//
//    indigox::Forcefield imp;
//
//    // Private wrapping functions
//    FFBond NewBondType(BondType t, int i, FFParam p) { return imp->NewBondType(t, i, p); }
//    FFAngle NewAngleType(AngleType t, int i, FFParam p) { return imp->NewAngleType(t, i, p); }
//    FFDihedral NewDihedralType(DihedralType t, int i, FFParam p) { return imp->NewDihedralType(t, i, p); }
//
//    // Public wrapping functions
//    TestForcefield() = delete;
//    TestForcefield(FFFamily f, std::string n) : imp(new IXForcefield(f,n)) { }
//    FFAtom NewAtomType(int i, std::string n) { return imp->NewAtomType(i,n); }
//    FFAtom GetAtomType(std::string n) const { return imp->GetAtomType(n); }
//    FFAtom GetAtomType(int i) const { return imp->GetAtomType(i); }
//    void ReserveAtomTypes(size_t s) { imp->ReserveAtomTypes(s); }
//    size_t NumAtomTypes() const { return imp->NumAtomTypes(); }
//
//    FFBond NewBondType(BondType t, int i, double a, double b) { return imp->NewBondType(t, i, a, b); }
//    FFBond GetBondType(BondType t, int i) const { return imp->GetBondType(t,i); }
//    FFBond GetBondType(int i) const { return imp->GetBondType(i); }
//    void LinkBondTypes(FFBond a, FFBond b) { imp->LinkBondTypes(a,b); }
//    void ReserveBondTypes(BondType t, size_t s) { imp->ReserveBondTypes(t, s); }
//    size_t NumBondTypes() const { return imp->NumBondTypes(); }
//    size_t NumBondTypes(BondType t) const { return imp->NumBondTypes(t); }
//
//    FFAngle NewAngleType(AngleType t, int i, double a, double b) { return imp->NewAngleType(t, i, a, b); }
//    FFAngle GetAngleType(AngleType t, int i) const { return imp->GetAngleType(t,i); }
//    FFAngle GetAngleType(int i) const { return imp->GetAngleType(i); }
//    void LinkAngleTypes(FFAngle a, FFAngle b) { imp->LinkAngleTypes(a,b); }
//    void ReserveAngleTypes(AngleType t, size_t s) { imp->ReserveAngleTypes(t, s); }
//    size_t NumAngleTypes() const { return imp->NumAngleTypes(); }
//    size_t NumAngleTypes(AngleType t) const { return imp->NumAngleTypes(t); }
//
//    FFDihedral NewDihedralType(DihedralType t, int i, double a, double b, double c) { return imp->NewDihedralType(t,i,a,b,c); }
//    FFDihedral NewDihedralType(DihedralType t, int i, double a, double b) { return imp->NewDihedralType(t,i,a,b); }
//    FFDihedral GetDihedralType(DihedralType t, int i) const { return imp->GetDihedralType(t,i); }
//    FFDihedral GetDihedralType(int i) const { return imp->GetDihedralType(i); }
//    void ReserveDihedralTypes(DihedralType t, size_t s) { imp->ReserveDihedralTypes(t,s); }
//    size_t NumDihedralTypes() const { return imp->NumDihedralTypes(); }
//    size_t NumDihedralTypes(DihedralType t) const { return imp->NumDihedralTypes(t); }
//
//    FFFamily GetFamily() const { return imp->GetFamily(); }
//    std::string GetName() const { return imp->GetName(); }
//
//    // internals acccess
//    FFFamily get_family() const { return imp->_family; }
//    std::string get_name() const { return imp->_name; }
//    AtmTypes& get_atms() { return imp->_atms; }
//    BndTypes& get_bnds() { return imp->_bnds; }
//    AngTypes& get_angs() { return imp->_angs; }
//    DhdTypes& get_dhds() { return imp->_dhds; }
  };
  
//  inline TestForcefield CreateGenericTestForcefield() {
//    return TestForcefield(FFFamily::GROMOS, "tester"); }
//
//  inline TestFFDihedral CreateGenericTestFFDihedral() {
//    return TestFFDihedral(DihedralType::Improper, 3, {0.95, 1.37},
//                          CreateGenericTestForcefield().imp); }
//  inline TestFFAngle CreateGenericTestFFAngle() {
//    return TestFFAngle(AngleType::Harmonic, 3, {0.56,0.0098},
//                       CreateGenericTestForcefield().imp); }
//  inline TestFFBond CreateGenericTestFFBond() {
//    return TestFFBond(BondType::Harmonic, 4, {0.567, 9.732},
//                      CreateGenericTestForcefield().imp); }
//  inline TestFFAtom CreateGenericTestFFAtom() {
//    return TestFFAtom(7, "GEN", CreateGenericTestForcefield().imp); }
  
  struct FFAtomTestFixture {
//    Forcefield ff;
//    TestFFAtom atm;
//    FFAtomTestFixture()
//    : ff(CreateGenericTestForcefield().imp), atm(7, "ATM", ff) { }
  };
  
  struct FFBondTestFixture {
//    Forcefield ff;
//    TestFFBond bnd, empty, harmonic, quartic;
//    FFBondTestFixture()
//    : ff(CreateGenericTestForcefield().imp), bnd(BondType::FENE, ff),
//    empty(BondType::Empty, 0, {}, ff),
//    harmonic(BondType::Harmonic, 3, {2.9227260e05, 1.09e-01}, ff),
//    quartic(BondType::Quartic, 3, {1.23e07, 1.09e-01}, ff) {
//      harmonic.set_link(quartic.imp);
//      quartic.set_link(harmonic.imp);
//    }
  };
  
  struct FFAngleTestFixture {
//    Forcefield ff;
//    TestFFAngle ang, empty, harmonic, cosineharmonic;
//    FFAngleTestFixture()
//    : ff(CreateGenericTestForcefield().imp), ang(AngleType::Quartic, ff),
//    empty(AngleType::Empty, 0, {}, ff),
//    harmonic(AngleType::Harmonic, 1, {1.1550101e-01, 9.0e+01}, ff),
//    cosineharmonic(AngleType::CosineHarmonic, 1, {3.8e+02, 9.0e+01}, ff) {
//      harmonic.set_link(cosineharmonic.imp);
//      cosineharmonic.set_link(harmonic.imp);
//    }
  };
  
  struct FFDihedralTestFixture {
//    Forcefield ff;
//    TestFFDihedral dhd, empty, proper, improper;
//    FFDihedralTestFixture()
//    : ff(CreateGenericTestForcefield().imp), dhd(DihedralType::Restricted, ff),
//    empty(DihedralType::Empty, 0, {}, ff),
//    proper(DihedralType::Proper, 1, {180.0, 2.67, 1}, ff),
//    improper(DihedralType::Improper, 2, {35.26439, 0.102}, ff) {  }
  };
  
  struct ForcefieldTestFixture {
    
  };
}

#endif
