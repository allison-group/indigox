//#include <algorithm>
#include <array>
#include <initializer_list>

#include <indigox/classes/forcefield.hpp>
#include <indigox/utils/quad.hpp>
#include <indigox/utils/triple.hpp>

#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/forcefield_test.hpp>

namespace indigox {
  
  test_suite_open("IXForcefield and types");
  
// ============================================================================
// == IXFFDihedral Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXFFDihedral::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("ff", _ff),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat));
  }
  
  template <class Archive>
  void IXFFDihedral::load_and_construct(Archive &archive,
                                        cereal::construct<IXFFDihedral> &construct,
                                        const uint32_t) {
    Type t;
    Forcefield ff;
    archive(INDIGOX_SERIAL_NVP("type", t),
            INDIGOX_SERIAL_NVP("ff", ff));
    
    construct(t, ff);
    archive(INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFDihedral);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFDihedral serilisation", T, ixffdhd_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::FFDihedralTestFixture fixture;
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(fixture.proper.imp, fixture.ff));
    }
    
    Forcefield loaded_ff;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(fixture.dhd.imp, loaded_ff));
    }
    
    check_eq(fixture.proper.get_id(), fixture.dhd.get_id());
    check_eq(fixture.proper.get_dat(), fixture.dhd.get_dat());
    check_eq(fixture.proper.get_mask(), fixture.dhd.get_mask());
    check_eq(fixture.proper.get_type(), fixture.dhd.get_type());
    check_eq(loaded_ff, fixture.dhd.GetForcefield());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffdhd_serial, ixserial<IXFFDihedral>);
  
// ============================================================================
// == IXFFDihedral Construction ===============================================
// ============================================================================
  
  IXFFDihedral::IXFFDihedral(DihedralType type, int id, FFParam parameters,
                             const Forcefield& ff)
  : IXFFDihedral(type, ff) {
    _id = id;
    if (parameters.size() != _mask.count())
      throw std::range_error("Incorrect item count");
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  IXFFDihedral::IXFFDihedral(DihedralType type, const Forcefield& ff)
  : _type(type), _ff(ff) {
    // Allow phase, force constant, multiplicty
    if (_type == DihedralType::Proper) _mask.from_uint32(7);
    // Allow force constant, ideal angle
    if (_type == DihedralType::Improper) _mask.from_uint32(10);
    _dat.fill(0);
  }
  
  test_case("IXFFDihedral construction") {
    Forcefield ff = test::CreateGenericTestForcefield().imp;
    using F = test::TestFFDihedral;
    using T = DihedralType;
    // Empty construction
    check_nothrow(F t(T::Empty, 4, {}, ff));
    check_nothrow(F t(T::Empty, ff));
    check_throws_as(F t(T::Empty, 4, {1.0}, ff), std::range_error);
    
    // Proper construction
    check_nothrow(F t(T::Proper, 4, {1.,2.,3}, ff));
    check_nothrow(F t(T::Proper, ff));
    check_throws_as(F t(T::Proper, 4, {1.,2.}, ff), std::range_error);
    check_throws_as(F t(T::Proper, 4, {1.,2.,3.,4.}, ff), std::range_error);
    
    // Improper construction
    check_nothrow(F t(T::Improper, 4, {1.,2.}, ff));
    check_nothrow(F t(T::Improper, ff));
    check_throws_as(F t(T::Improper, 4, {1.}, ff), std::range_error);
    check_throws_as(F t(T::Improper, 4, {1.,2.,3.}, ff), std::range_error);
  }
  
  test_case_fixture(test::FFDihedralTestFixture, "IXFFDihedral getting checks") {
    // Empty
    check_throws_as(empty.GetPhaseShift(), std::runtime_error);
    check_throws_as(empty.GetForceConstant(), std::runtime_error);
    check_throws_as(empty.GetMultiplicity(), std::runtime_error);
    check_throws_as(empty.GetIdealAngle(), std::runtime_error);
    check_eq(DihedralType::Empty, empty.GetType());
    check_eq(0, empty.GetID());
    check_eq(ff, empty.GetForcefield());
    
    // Proper
    check_eq(approximately(180.0), proper.GetPhaseShift());
    check_eq(approximately(2.67), proper.GetForceConstant());
    check_eq(1, proper.GetMultiplicity());
    check_throws_as(proper.GetIdealAngle(), std::runtime_error);
    check_eq(DihedralType::Proper, proper.GetType());
    check_eq(1, proper.GetID());
    check_eq(ff, proper.GetForcefield());
    
    // Improper
    check_throws_as(improper.GetPhaseShift(), std::runtime_error);
    check_eq(approximately(0.102), improper.GetForceConstant());
    check_throws_as(improper.GetMultiplicity(), std::runtime_error);
    check_eq(approximately(35.26439), improper.GetIdealAngle());
    check_eq(DihedralType::Improper, improper.GetType());
    check_eq(2, improper.GetID());
    check_eq(ff, improper.GetForcefield());
  }
  
  test_case_fixture(test::FFDihedralTestFixture, "IXFFDihedral internals checks") {
    using AllowEnum = indigox::test::TestFFDihedral::AllowEnum;
    // Empty
    check_eq(DihedralType::Empty, empty.get_type());
    check_eq(0, empty.get_id());
    check_eq(approximately(0.0), empty.get_dat()[0]);
    check_eq(approximately(0.0), empty.get_dat()[1]);
    check_eq(approximately(0.0), empty.get_dat()[2]);
    check_false(empty.get_mask()[AllowEnum::Allow_IdealAngle]);
    check_false(empty.get_mask()[AllowEnum::Allow_PhaseShift]);
    check_false(empty.get_mask()[AllowEnum::Allow_Multiplicity]);
    check_false(empty.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(ff, empty.get_ff().lock());
    
    // Proper
    check_eq(DihedralType::Proper, proper.get_type());
    check_eq(1, proper.get_id());
    check_eq(approximately(180.0), proper.get_dat()[0]);
    check_eq(approximately(2.67), proper.get_dat()[1]);
    check_eq(approximately(1.0), proper.get_dat()[2]);
    check_false(proper.get_mask()[AllowEnum::Allow_IdealAngle]);
    check(proper.get_mask()[AllowEnum::Allow_PhaseShift]);
    check(proper.get_mask()[AllowEnum::Allow_Multiplicity]);
    check(proper.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(ff, proper.get_ff().lock());
    
    // Improper
    check_eq(DihedralType::Improper, improper.get_type());
    check_eq(2, improper.get_id());
    check_eq(approximately(35.26439), improper.get_dat()[0]);
    check_eq(approximately(0.102), improper.get_dat()[1]);
    check_eq(approximately(0.0), improper.get_dat()[2]);
    check(improper.get_mask()[AllowEnum::Allow_IdealAngle]);
    check_false(improper.get_mask()[AllowEnum::Allow_PhaseShift]);
    check_false(improper.get_mask()[AllowEnum::Allow_Multiplicity]);
    check(improper.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(ff, improper.get_ff().lock());
  }
  
// ============================================================================
// == IXFFAngle Serialisation =================================================
// ============================================================================
  
  template <class Archive>
  void IXFFAngle::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("ff", _ff),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("link", _link));
  }
  
  template <class Archive>
  void IXFFAngle::load_and_construct(Archive &archive,
                                            cereal::construct<IXFFAngle> &construct,
                                            const uint32_t) {
    Type t;
    Forcefield ff;
    archive(INDIGOX_SERIAL_NVP("type", t),
            INDIGOX_SERIAL_NVP("ff", ff));
    
    construct(t, ff);
    archive(INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("link", construct->_link));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFAngle);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAngle serilisation", T, ixffang_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::FFAngleTestFixture fixture;
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(fixture.harmonic.imp, fixture.ff));
    }
    
    Forcefield loaded_ff;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(fixture.ang.imp, loaded_ff));
    }
    
    check_eq(fixture.harmonic.get_id(), fixture.ang.get_id());
    check_eq(fixture.harmonic.get_dat(), fixture.ang.get_dat());
    check_eq(fixture.harmonic.get_mask(), fixture.ang.get_mask());
    check_eq(fixture.harmonic.get_type(), fixture.ang.get_type());
    check_eq(loaded_ff, fixture.ang.GetForcefield());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffang_serial, ixserial<IXFFAngle>);
  
// ============================================================================
// == IXFFAngle Construction ==============================================
// ============================================================================
  
  IXFFAngle::IXFFAngle(AngleType type, int id, FFParam parameters,
                       const Forcefield& ff)
  : IXFFAngle(type, ff) {
    _id = id;
    if (parameters.size() != _mask.count())
      throw std::range_error("Incorrect item count");
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  IXFFAngle::IXFFAngle(AngleType type, const Forcefield& ff)
  : _type(type), _ff(ff) {
    // Allow force constant and ideal angle
    if (_type == AngleType::Harmonic) _mask.from_uint32(3);
    // Allow force constant and ideal angle
    if (_type == AngleType::CosineHarmonic) _mask.from_uint32(3);
    _dat.fill(0);
  }
  
  test_case("IXFFAngle construction") {
    Forcefield ff = test::CreateGenericTestForcefield().imp;
    using F = test::TestFFAngle;
    using T = AngleType;
    // Empty construction
    check_nothrow(F t(T::Empty, 4, {}, ff));
    check_nothrow(F t(T::Empty, ff));
    check_throws_as(F t(T::Empty, 4, {1.}, ff), std::range_error);
    
    // Harmonic construction
    check_nothrow(F t(T::Harmonic, 4, {1.,2.}, ff));
    check_nothrow(F t(T::Harmonic, ff));
    check_throws_as(F t(T::Harmonic, 4, {1.}, ff), std::range_error);
    check_throws_as(F t(T::Empty, 4, {1.,2.,3.}, ff), std::range_error);
    
    // CosineHarmonic construction
    check_nothrow(F t(T::CosineHarmonic, 4, {1.,2.}, ff));
    check_nothrow(F t(T::CosineHarmonic, ff));
    check_throws_as(F t(T::CosineHarmonic, 4, {1.}, ff), std::range_error);
    check_throws_as(F t(T::CosineHarmonic, 4, {1.,2.,3.}, ff), std::range_error);
  }
  
  test_case_fixture(test::FFAngleTestFixture, "IXFFAngle getting checks") {
    // Empty
    check_throws_as(empty.GetIdealAngle(), std::runtime_error);
    check_throws_as(empty.GetForceConstant(), std::runtime_error);
    check_eq(AngleType::Empty, empty.GetType());
    check_eq(0, empty.GetID());
    check_eq(ff, empty.GetForcefield());
    
    // Harmonic
    check_eq(approximately(90.), harmonic.GetIdealAngle());
    check_eq(approximately(0.11550101), harmonic.GetForceConstant());
    check_eq(AngleType::Harmonic, harmonic.GetType());
    check_eq(1, harmonic.GetID());
    check_eq(ff, harmonic.GetForcefield());
    check_eq(cosineharmonic.imp, harmonic.GetLinkedType());
    
    // CosineHarmonic
    check_eq(approximately(90.), cosineharmonic.GetIdealAngle());
    check_eq(approximately(380.), cosineharmonic.GetForceConstant());
    check_eq(AngleType::CosineHarmonic, cosineharmonic.GetType());
    check_eq(1, cosineharmonic.GetID());
    check_eq(ff, cosineharmonic.GetForcefield());
    check_eq(harmonic.imp, cosineharmonic.GetLinkedType());
  }
  
  test_case_fixture(test::FFAngleTestFixture, "IXFFAngle internals checks") {
    using AllowEnum = indigox::test::TestFFAngle::AllowEnum;
    // Empty
    check_eq(AngleType::Empty, empty.get_type());
    check_eq(0, empty.get_id());
    check_eq(approximately(0.0), empty.get_dat()[0]);
    check_eq(approximately(0.0), empty.get_dat()[1]);
    check_false(empty.get_mask()[AllowEnum::Allow_IdealAngle]);
    check_false(empty.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(FFAngle(), empty.get_link());
    check_eq(ff, empty.get_ff().lock());
    
    // Harmonic
    check_eq(AngleType::Harmonic, harmonic.get_type());
    check_eq(1, harmonic.get_id());
    check_eq(approximately(0.11550101), harmonic.get_dat()[0]);
    check_eq(approximately(90.), harmonic.get_dat()[1]);
    check(harmonic.get_mask()[AllowEnum::Allow_IdealAngle]);
    check(harmonic.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(cosineharmonic.imp, harmonic.get_link());
    check_eq(ff, harmonic.get_ff().lock());
    
    // CosineHarmonic
    check_eq(AngleType::CosineHarmonic, cosineharmonic.get_type());
    check_eq(1, cosineharmonic.get_id());
    check_eq(approximately(380.), cosineharmonic.get_dat()[0]);
    check_eq(approximately(90.), cosineharmonic.get_dat()[1]);
    check(cosineharmonic.get_mask()[AllowEnum::Allow_IdealAngle]);
    check(cosineharmonic.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(harmonic.imp, cosineharmonic.get_link());
    check_eq(ff, cosineharmonic.get_ff().lock());
  }

// ============================================================================
// == IXFFBond Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXFFBond::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("ff", _ff),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("link", _link));
  }
  
  template <class Archive>
  void IXFFBond::load_and_construct(Archive &archive,
                                        cereal::construct<IXFFBond> &construct,
                                        const uint32_t) {
    Type t;
    Forcefield ff;
    archive(INDIGOX_SERIAL_NVP("type", t),
            INDIGOX_SERIAL_NVP("ff", ff));
    
    construct(t, ff);
    archive(INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("link", construct->_link));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFBond);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFBond serilisation", T, ixffbnd_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::FFBondTestFixture fixture;
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(fixture.harmonic.imp, fixture.ff));
    }
    
    Forcefield loaded_ff;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(fixture.bnd.imp, loaded_ff));
    }
    
    check_eq(fixture.harmonic.get_id(), fixture.bnd.get_id());
    check_eq(fixture.harmonic.get_dat(), fixture.bnd.get_dat());
    check_eq(fixture.harmonic.get_mask(), fixture.bnd.get_mask());
    check_eq(fixture.harmonic.get_type(), fixture.bnd.get_type());
    check_eq(loaded_ff, fixture.bnd.GetForcefield());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffbnd_serial, ixserial<IXFFBond>);
  
// ============================================================================
// == IXFFBond Construction ===============================================
// ============================================================================
  
  IXFFBond::IXFFBond(BondType type, int id, FFParam parameters,
                     const Forcefield& ff)
  : IXFFBond(type, ff) {
    _id = id;
    if (parameters.size() != _mask.count())
      throw std::range_error("Incorrect item count");
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  IXFFBond::IXFFBond(BondType type, const Forcefield& ff)
  : _type(type), _ff(ff) {
    // Allow force constant and ideal length
    if (_type == BondType::Harmonic) _mask.from_uint32(3);
    // Allow force constant and ideal length
    if (_type == BondType::Quartic) _mask.from_uint32(3);
    _dat.fill(0);
  }
  
  test_case("IXFFBond construction") {
    Forcefield ff = test::CreateGenericTestForcefield().imp;
    using F = test::TestFFBond;
    using T = BondType;
    // Empty construction
    check_nothrow(F t(T::Empty, 4, {}, ff));
    check_nothrow(F t(T::Empty, ff));
    check_throws_as(F t(T::Empty, 4, {1.}, ff), std::range_error);
    
    // Harmonic construction
    check_nothrow(F t(T::Harmonic, 4, {1.,2.}, ff));
    check_nothrow(F t(T::Harmonic, ff));
    check_throws_as(F t(T::Harmonic, 4, {1.}, ff), std::range_error);
    check_throws_as(F t(T::Empty, 4, {1.,2.,3.}, ff), std::range_error);
    
    // Quartic construction
    check_nothrow(F t(T::Quartic, 4, {1.,2.}, ff));
    check_nothrow(F t(T::Quartic, ff));
    check_throws_as(F t(T::Quartic, 4, {1.}, ff), std::range_error);
    check_throws_as(F t(T::Quartic, 4, {1.,2.,3.}, ff), std::range_error);
  }
  
  test_case_fixture(test::FFBondTestFixture, "IXFFBond getting checks") {
    // Empty
    check_throws_as(empty.GetIdealLength(), std::runtime_error);
    check_throws_as(empty.GetForceConstant(), std::runtime_error);
    check_eq(BondType::Empty, empty.GetType());
    check_eq(0, empty.GetID());
    check_eq(ff, empty.GetForcefield());
    
    // Harmonic
    check_eq(approximately(0.109), harmonic.GetIdealLength());
    check_eq(approximately(292272.6), harmonic.GetForceConstant());
    check_eq(BondType::Harmonic, harmonic.GetType());
    check_eq(3, harmonic.GetID());
    check_eq(ff, harmonic.GetForcefield());
    check_eq(quartic.imp, harmonic.GetLinkedType());
    
    // Quartic
    check_eq(approximately(0.109), quartic.GetIdealLength());
    check_eq(approximately(12300000.), quartic.GetForceConstant());
    check_eq(BondType::Quartic, quartic.GetType());
    check_eq(3, quartic.GetID());
    check_eq(ff, quartic.GetForcefield());
    check_eq(harmonic.imp, quartic.GetLinkedType());
  }
  
  test_case_fixture(test::FFBondTestFixture, "IXFFBond internals checks") {
    using AllowEnum = indigox::test::TestFFBond::AllowEnum;
    // Empty
    check_eq(BondType::Empty, empty.get_type());
    check_eq(0, empty.get_id());
    check_eq(approximately(0.0), empty.get_dat()[0]);
    check_eq(approximately(0.0), empty.get_dat()[1]);
    check_false(empty.get_mask()[AllowEnum::Allow_IdealLength]);
    check_false(empty.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(FFBond(), empty.get_link());
    check_eq(ff, empty.get_ff().lock());
    
    // Harmonic
    check_eq(BondType::Harmonic, harmonic.get_type());
    check_eq(3, harmonic.get_id());
    check_eq(approximately(292272.6), harmonic.get_dat()[0]);
    check_eq(approximately(0.109), harmonic.get_dat()[1]);
    check(harmonic.get_mask()[AllowEnum::Allow_IdealLength]);
    check(harmonic.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(quartic.imp, harmonic.get_link());
    check_eq(ff, harmonic.get_ff().lock());
    
    // Quartic
    check_eq(BondType::Quartic, quartic.get_type());
    check_eq(3, quartic.get_id());
    check_eq(approximately(12300000.), quartic.get_dat()[0]);
    check_eq(approximately(0.109), quartic.get_dat()[1]);
    check(quartic.get_mask()[AllowEnum::Allow_IdealLength]);
    check(quartic.get_mask()[AllowEnum::Allow_ForceConstant]);
    check_eq(harmonic.imp, quartic.get_link());
    check_eq(ff, quartic.get_ff().lock());
  }
  
// ============================================================================
// == IXFFAtom Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXFFAtom::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("name", _name),
            INDIGOX_SERIAL_NVP("ff", _ff));
  }
  
  template <class Archive>
  void IXFFAtom::load_and_construct(Archive &archive,
                                    cereal::construct<IXFFAtom> &construct,
                                    const uint32_t) {
    int id;
    std::string name;
    Forcefield ff;
    archive(INDIGOX_SERIAL_NVP("id", id),
            INDIGOX_SERIAL_NVP("name", name),
            INDIGOX_SERIAL_NVP("ff", ff));
    construct(id, name, ff);
  }
  INDIGOX_SERIALISE_SPLIT(IXFFAtom);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAtom serialisation", T, ixffatm_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::FFAtomTestFixture fixture;
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(fixture.atm.imp, fixture.ff));
    }
    
    test::TestFFAtom loaded(0, "", fixture.ff);
    Forcefield loaded_ff;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded.imp, loaded_ff));
    }
    
    check_eq(fixture.atm.get_id(), loaded.get_id());
    check_eq(fixture.atm.get_name(), loaded.get_name());
    check_eq(loaded_ff, loaded.get_ff().lock());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffatm_serial, ixserial<IXFFAtom>);
  
// ============================================================================
// == IXFFAtom Construction ===============================================
// ============================================================================
  
  IXFFAtom::IXFFAtom(int id, std::string name, const Forcefield& ff)
  : _id(id), _name(name), _ff(ff) { }
  
  test_case("IXFFAtom construction") {
    Forcefield ff = test::CreateGenericTestForcefield().imp;
    using F = test::TestFFAtom;
    check_nothrow(F t(23, "TEST", ff));
  }
  
  test_case_fixture(test::FFAtomTestFixture, "IXFFAtom getting checks") {
    check_eq(7, atm.GetID());
    check_eq("ATM", atm.GetName());
    check_eq(ff, atm.GetForcefield());
  }
  
  test_case_fixture(test::FFAtomTestFixture, "IXFFAtom internals checks") {
    check_eq(7, atm.get_id());
    check_eq("ATM", atm.get_name());
    check_eq(ff, atm.get_ff().lock());
  }
  
// ============================================================================
// == IXForcefield Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXForcefield::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("family", _family),
            INDIGOX_SERIAL_NVP("name", _name),
            INDIGOX_SERIAL_NVP("atom_types", _atms),
            INDIGOX_SERIAL_NVP("bond_types", _bnds),
            INDIGOX_SERIAL_NVP("angle_types", _angs),
            INDIGOX_SERIAL_NVP("dihedral_types", _dhds));
  }
  
  template <class Archive>
  void IXForcefield::load_and_construct(Archive &archive,
                                        cereal::construct<IXForcefield> &construct,
                                        const uint32_t) {
    construct(FFFamily::Empty, "");
    archive(INDIGOX_SERIAL_NVP("family", construct->_family),
            INDIGOX_SERIAL_NVP("name", construct->_name),
            INDIGOX_SERIAL_NVP("atom_types", construct->_atms),
            INDIGOX_SERIAL_NVP("bond_types", construct->_bnds),
            INDIGOX_SERIAL_NVP("angle_types", construct->_angs),
            INDIGOX_SERIAL_NVP("dihedral_types", construct->_dhds));
  }
  
// ============================================================================
// == IXForcefield Construction ===============================================
// ============================================================================
  
  IXForcefield::IXForcefield(FFFamily family, std::string name)
  : _family(family), _name(name) {
    
    if (family == FFFamily::GROMOS) {
      _bnds.emplace(BondType::Harmonic, BondTypes::mapped_type());
      _bnds.emplace(BondType::Quartic, BondTypes::mapped_type());
      _angs.emplace(AngleType::Harmonic, AngleTypes::mapped_type());
      _angs.emplace(AngleType::CosineHarmonic, AngleTypes::mapped_type());
      _dhds.emplace(DihedralType::Proper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralType::Improper, DihedralTypes::mapped_type());
    }
    if (family == FFFamily::Other) {
      _bnds.emplace(BondType::Harmonic, BondTypes::mapped_type());
      _bnds.emplace(BondType::Quartic, BondTypes::mapped_type());
      _bnds.emplace(BondType::Morse, BondTypes::mapped_type());
      _bnds.emplace(BondType::Cubic, BondTypes::mapped_type());
      _bnds.emplace(BondType::FENE, BondTypes::mapped_type());
      _angs.emplace(AngleType::Harmonic, AngleTypes::mapped_type());
      _angs.emplace(AngleType::CosineHarmonic, AngleTypes::mapped_type());
      _angs.emplace(AngleType::UreyBradley, AngleTypes::mapped_type());
      _angs.emplace(AngleType::Quartic, AngleTypes::mapped_type());
      _dhds.emplace(DihedralType::Proper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralType::Improper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralType::RyckaertBellemans, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralType::PeriodicImproper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralType::Fourier, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralType::Restricted, DihedralTypes::mapped_type());
    }
  }
  
// ============================================================================
// == IXForcefield Adding types ===============================================
// ============================================================================
  
  FFAtom IXForcefield::NewAtomType(int id, std::string name) {
    if (GetAtomType(id) || GetAtomType(name)) return FFAtom();
    _atms.emplace_back(std::make_shared<IXFFAtom>(id, name, shared_from_this()));
    return _atms.back();
  }
  
  FFBond IXForcefield::NewBondType(BondType type, int id, FFParam param) {
    if (_bnds.find(type) == _bnds.end()) return FFBond();
    if (GetBondType(type, id)) return FFBond();
    _bnds[type].emplace_back(std::make_shared<IXFFBond>(type, id, param,
                                                        shared_from_this()));
    return _bnds[type].back();
  }
  
  FFAngle IXForcefield::NewAngleType(AngleType type, int id, FFParam param) {
    if (_angs.find(type) == _angs.end()) return FFAngle();
    if (GetAngleType(type, id)) return FFAngle();
    _angs[type].emplace_back(std::make_shared<IXFFAngle>(type, id, param,
                                                         shared_from_this()));
    return _angs[type].back();
  }
  
  FFDihedral IXForcefield::NewDihedralType(DihedralType type, int id, FFParam param) {
    if (_dhds.find(type) == _dhds.end()) return FFDihedral();
    if (GetDihedralType(type, id)) return FFDihedral();
    _dhds[type].emplace_back(std::make_shared<IXFFDihedral>(type, id, param,
                                                            shared_from_this()));
    return _dhds[type].back();
  }
  
// ============================================================================
// == IXForcefield Atom type handlings ========================================
// ============================================================================
  
  FFAtom IXForcefield::GetAtomType(std::string name) const {
    auto pred = [&name](const FFAtom& t) { return t->GetName() == name; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? FFAtom() : *pos;
  }
  
  FFAtom IXForcefield::GetAtomType(int id) const {
    auto pred = [&id](const FFAtom& t) { return t->GetID() == id; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? FFAtom() : *pos;
  }

// ============================================================================
// == IXForcefield Bond type handlings ========================================
// ============================================================================
  
  FFBond IXForcefield::GetBondType(BondType type, int id) const {
    if (_bnds.find(type) == _bnds.end()) return FFBond();
    auto pred = [&id](const FFBond& t) { return t->GetID() == id; };
    auto pos = std::find_if(_bnds.at(type).begin(), _bnds.at(type).end(), pred);
    return (pos == _bnds.at(type).end()) ? FFBond() : *pos;
  }
  
  FFBond IXForcefield::GetBondType(int id) const {
    for (auto& types: _bnds) {
      FFBond found = GetBondType(types.first, id);
      if (found) return found;
    }
    return FFBond();
  }
  
  void IXForcefield::LinkBondTypes(FFBond a, FFBond b) {
    if (a->GetLinkedType()) a->GetLinkedType()->_link.reset();
    if (b->GetLinkedType()) b->GetLinkedType()->_link.reset();
    a->_link = b;
    b->_link = a;
  }
  
// ============================================================================
// == IXForcefield Angle type handlings =======================================
// ============================================================================
  
  FFAngle IXForcefield::GetAngleType(AngleType type, int id) const {
    if (_angs.find(type) == _angs.end()) return FFAngle();
    auto pred = [&id](const FFAngle& t) { return t->GetID() == id; };
    auto pos = std::find_if(_angs.at(type).begin(), _angs.at(type).end(), pred);
    return (pos == _angs.at(type).end()) ? FFAngle() : *pos;
  }
  
  FFAngle IXForcefield::GetAngleType(int id) const {
    for (auto& types: _angs) {
      FFAngle found = GetAngleType(types.first, id);
      if (found) return found;
    }
    return FFAngle();
  }
  
  void IXForcefield::LinkAngleTypes(FFAngle a, FFAngle b) {
    if (a->GetLinkedType()) a->GetLinkedType()->_link.reset();
    if (b->GetLinkedType()) b->GetLinkedType()->_link.reset();
    a->_link = b;
    b->_link = a;
  }
  
// ============================================================================
// == IXForcefield Dihedral type handlings ====================================
// ============================================================================
  
  FFDihedral IXForcefield::GetDihedralType(DihedralType type, int id) const {
    if (_dhds.find(type) == _dhds.end()) return FFDihedral();
    auto pred = [&id](const FFDihedral& t) { return t->GetID() == id; };
    auto pos = std::find_if(_dhds.at(type).begin(), _dhds.at(type).end(), pred);
    return (pos == _dhds.at(type).end()) ? FFDihedral() : *pos;
  }
  
  FFDihedral IXForcefield::GetDihedralType(int id) const {
    for (auto& types: _dhds) {
      FFDihedral found = GetDihedralType(types.first, id);
      if (found) return found;
    }
    return FFDihedral();
  }
  
  
// ============================================================================
// == Hardcoded forcefield generations ========================================
// ============================================================================
  Forcefield GenerateGROMOS54A7() {
    Forcefield ff = std::make_shared<IXForcefield>(FFFamily::GROMOS, "54A7");
    // Add atom types
    std::vector<std::pair<int, std::string>> atom_dat = {
      {2, "OM"}, {1, "O"},
      {4, "OE"}, {39, "CChl"}, {32, "F"}, {52, "OUrea"}, {15, "CH2"}, {20, "HC"},
      {37, "NA+"}, {33, "CL"}, {44, "ODmso"}, {40, "CLChl"}, {45, "CCl4"},
      {31, "AR"}, {38, "CL-"}, {35, "CMet"}, {21, "H"}, {27, "ZN2+"}, {28, "MG2+"},
      {46, "CLCl4"}, {36, "OMet"}, {54, "CH3p"}, {51, "CUrea"}, {30, "P,SI"},
      {7, "NT"}, {26, "FE"}, {43, "CDmso"}, {17, "CH4"}, {29, "CA2+"}, {22, "DUM"},
      {42, "SDmso"}, {48, "CTFE"}, {5, "OW"}, {34, "BR"}, {49, "CHTFE"}, {14, "CH1"},
      {41, "HChl"}, {53, "NUrea"}, {23, "S"}, {13, "CH0"}, {19, "CR1"}, {11, "NE"},
      {16, "CH3"}, {8, "NL"}, {12, "C"}, {10, "NZ"}, {6, "N"}, {24, "CU1+"},
      {50, "OTFE"}, {25, "CU2+"}, {3, "OA"}, {47, "FTFE"}, {18, "CH2r"}, {9, "NR"}};
    for (auto& id_name : atom_dat) ff->NewAtomType(id_name.first, id_name.second);
    
    // Add bond types
    std::vector<stdx::quad<int, double, double, double>> bnd_dat = {
      {34, 0.198, 50181.12, 640000.0}, {48, 0.290283, 502214.75, 2980000.0},
      {22, 0.148, 251019.84, 5730000.0}, {49, 0.279388, 373115.59, 2390000.0},
      {50, 0.291189, 371384.73, 2190000.0}, {10, 0.133, 417460.4, 11800000.0},
      {23, 0.148, 334693.12, 7640000.0}, {4, 0.112, 928256.0, 37000000.0},
      {14, 0.138, 418968.0, 11000000.0}, {31, 0.178, 376405.92, 5940000.0},
      {27, 0.153, 334748.7, 7150000.0}, {35, 0.2, 50240.0, 628000.0},
      {6, 0.125, 418750.0, 13400000.0}, {5, 0.123, 502282.8, 16600000.0},
      {7, 0.132, 418176.0, 12000000.0}, {12, 0.134, 420170.4, 11700000.0},
      {39, 0.11, 292820.0, 12100000.0}, {19, 0.143, 376670.58, 9210000.0},
      {28, 0.161, 250915.28, 4840000.0}, {21, 0.147, 376428.78, 8710000.0},
      {37, 0.221, 52748.28, 540000.0}, {46, 0.163299, 464531.53, 8710000.0},
      {45, 0.135, 375435.0, 10300000.0}, {16, 0.139, 417333.6, 10800000.0},
      {52, 0.287407, 502224.92, 3040000.0}, {41, 0.153, 376416.72, 8040000.0},
      {13, 0.136, 377318.4, 10200000.0}, {24, 0.148, 376748.8, 8600000.0},
      {25, 0.15, 376650.0, 8370000.0}, {17, 0.14, 334768.0, 8540000.0},
      {2, 0.1, 374000.0, 18700000.0}, {33, 0.187, 251077.42, 3590000.0},
      {36, 0.204, 418656.96, 5030000.0}, {15, 0.139, 334639.72, 8660000.0},
      {20, 0.1435, 251225.45, 6100000.0}, {43, 0.176, 501811.2, 8100000.0},
      {40, 0.1758, 501907.59, 8120000.0}, {32, 0.183, 376416.36, 5620000.0},
      {1, 0.1, 314000.0, 15700000.0}, {3, 0.109, 292272.6, 12300000.0},
      {30, 0.178, 172360.96, 2720000.0}, {11, 0.134, 377076.0, 10500000.0},
      {9, 0.133, 375006.8, 10600000.0}, {47, 0.233839, 293088.43, 2680000.0},
      {26, 0.152, 250909.44, 5430000.0}, {8, 0.133, 313802.86, 8870000.0},
      {42, 0.193799, 371824.72, 4950000.0}, {18, 0.143, 334545.64, 8180000.0},
      {29, 0.163, 250811.36, 4720000.0}, {38, 0.1, 464000.0, 23200000.0},
      {44, 0.1265, 419258.95, 13100000.0}, {51, 0.2077, 342525.96, 3970000.0}};
    for (auto& dat : bnd_dat)
      ff->LinkBondTypes(ff->NewBondType(BondType::Harmonic, dat.first, dat.third, dat.second),
                        ff->NewBondType(BondType::Quartic, dat.first, dat.fourth, dat.second));
    
    // Add angle types
    std::vector<stdx::quad<int, double, double, double>> ang_dat = {
      {1, 90.0, 0.11550101, 380.0}, {3, 96.0, 0.12177061, 405.0},
      {14, 109.6, 0.12142334, 450.0}, {9, 109.5, 0.08638614, 320.0},
      {40, 155.0, 0.12112698, 2215.0}, {15, 111.0, 0.14048747, 530.0},
      {29, 120.0, 0.17801113, 780.0}, {34, 125.0, 0.076490216, 375.0},
      {24, 120.0, 0.10147593, 445.0}, {37, 126.0, 0.12744672, 640.0},
      {45, 97.4, 0.14024534, 469.0}, {28, 120.0, 0.15288018, 670.0},
      {20, 116.0, 0.11421859, 465.0}, {12, 109.5, 0.12157397, 450.0},
      {27, 120.0, 0.12774922, 560.0}, {36, 126.0, 0.11448735, 575.0},
      {43, 107.57, 0.13376523, 484.0}, {46, 106.75, 0.14026005, 503.0},
      {51, 110.3, 0.14017954, 524.0}, {5, 103.0, 0.12122177, 420.0},
      {8, 109.5, 0.076912479, 285.0}, {47, 108.53, 0.12108416, 443.0},
      {19, 115.0, 0.1524165, 610.0}, {41, 180.0, 0.072640156, 91350.0},
      {16, 113.0, 0.14045138, 545.0}, {4, 100.0, 0.14008261, 475.0},
      {2, 90.0, 0.12768574, 420.0}, {50, 109.5, 0.12103261, 448.0},
      {18, 115.0, 0.11488482, 460.0}, {17, 115.0, 0.01229657, 50.0},
      {32, 123.0, 0.088743846, 415.0}, {44, 111.3, 0.16689058, 632.0},
      {39, 132.0, 0.12775497, 760.0}, {13, 109.5, 0.14052124, 520.0},
      {31, 122.0, 0.15317431, 700.0}, {53, 117.2, 0.15305438, 636.0},
      {42, 109.5, 0.11724316, 434.0}, {21, 116.0, 0.15236094, 620.0},
      {30, 121.0, 0.15312732, 685.0}, {35, 125.0, 0.1531408, 750.0},
      {6, 104.0, 0.14028506, 490.0}, {48, 109.5, 0.1670474, 618.0},
      {52, 111.4, 0.14025677, 532.0}, {54, 121.4, 0.1529482, 690.0},
      {7, 108.0, 0.12788754, 465.0}, {22, 117.0, 0.15336019, 635.0},
      {11, 109.5, 0.11480708, 425.0}, {26, 120.0, 0.12089532, 530.0},
      {10, 109.5, 0.10262668, 380.0}, {33, 124.0, 0.15266919, 730.0},
      {25, 120.0, 0.11518373, 505.0}, {38, 126.0, 0.15336544, 770.0},
      {49, 107.6, 0.14008648, 507.0}, {23, 120.0, 0.088910434, 390.0}};
    for (auto& dat : ang_dat)
      ff->LinkAngleTypes(ff->NewAngleType(AngleType::Harmonic, dat.first,
                                          dat.third, dat.second),
                         ff->NewAngleType(AngleType::CosineHarmonic, dat.first,
                                          dat.fourth, dat.second));
    
    // Add improper types
    std::vector<stdx::triple<int, double, double>> imp_dat = {
      {4, 0.051, 180.0}, {2, 0.102, 35.26439}, {5, 0.102, -35.26439},
      {3, 0.204, 0.0}, {1, 0.051, 0.0}};
    for (auto& dat : imp_dat)
      ff->NewDihedralType(DihedralType::Improper, dat.first, dat.second, dat.third);
    
    // Add proper types
    std::vector<stdx::quad<int, double, double, int>> prp_dat = {
      {19, 3.14, 0.0, 2}, {42, 3.5, 180.0, 2}, {17, 0.418, 0.0, 2},
      {2, 3.41, 180.0, 1}, {7, 2.79, 0.0, 1}, {24, 1.3, 0.0, 3},
      {37, 9.5, 0.0, 3}, {45, 0.4, 0.0, 6}, {15, 41.8, 180.0, 2},
      {16, 0.0, 0.0, 2}, {38, 0.0, 0.0, 4}, {41, 3.77, 0.0, 6},
      {44, 0.7, 180.0, 6}, {22, 1.05, 0.0, 3}, {6, 9.45, 180.0, 1},
      {10, 5.86, 180.0, 2}, {28, 3.65, 0.0, 3}, {13, 24.0, 180.0, 2},
      {1, 2.67, 180.0, 1}, {35, 7.69, 0.0, 3}, {11, 7.11, 180.0, 2},
      {30, 3.9, 0.0, 3}, {3, 4.97, 180.0, 1}, {20, 5.09, 0.0, 2},
      {23, 1.26, 0.0, 3}, {32, 4.69, 0.0, 3}, {26, 2.93, 0.0, 3},
      {25, 2.53, 0.0, 3}, {9, 1.53, 180.0, 2}, {12, 16.7, 180.0, 2},
      {27, 3.19, 0.0, 3}, {14, 33.5, 180.0, 2}, {8, 5.35, 0.0, 1},
      {33, 5.44, 0.0, 3}, {34, 5.92, 0.0, 3}, {5, 9.35, 180.0, 1},
      {29, 3.77, 0.0, 3}, {39, 1.0, 180.0, 6}, {4, 5.86, 180.0, 1},
      {31, 4.18, 0.0, 3}, {43, 2.8, 0.0, 3}, {36, 8.62, 0.0, 3},
      {40, 1.0, 0.0, 6}, {21, 16.7, 0.0, 2}, {18, 2.09, 0.0, 2}};
    for (auto& dat : prp_dat)
      ff->NewDihedralType(DihedralType::Proper, dat.first, dat.second,
                          dat.third, dat.fourth);
    
    return ff;
  }
  
  test_suite_close();
}
