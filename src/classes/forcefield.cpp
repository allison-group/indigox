//#include <algorithm>
#include <array>
#include <initializer_list>

#include <indigox/classes/forcefield.hpp>
#include <indigox/utils/numerics.hpp>

#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/forcefield_test.hpp>

namespace indigox {
  
  test_suite_open("IXForcefield and types");
  
// ============================================================================
// == IXFFDihedral Serialisation ==========================================
// ============================================================================
  
  template <class Archive>
  void IXFFDihedral::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("mask", _mask));
  }
  
  template <class Archive>
  void IXFFDihedral::load_and_construct(Archive &archive,
                                            cereal::construct<IXFFDihedral> &construct,
                                            const uint32_t) {
    construct(DihedralType::Empty, 0, std::initializer_list<float_>());
    archive(INDIGOX_SERIAL_NVP("type", construct->_type),
            INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("mask", construct->_mask));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFDihedral);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFDihedral serilisation", T, ixffdhd_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestFFDihedral saved(DihedralType::Proper, 72, {0.1,0.2,0.3});
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestFFDihedral loaded(DihedralType::Empty, 0, {});
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded.imp));
    }
    
    check_eq(saved.get_id(), loaded.get_id());
    check_eq(saved.get_dat(), loaded.get_dat());
    check_eq(saved.get_mask(), loaded.get_mask());
    check_eq(saved.get_type(), loaded.get_type());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffdhd_serial, ixserial<IXFFDihedral>);
  
// ============================================================================
// == IXFFDihedral Construction ===========================================
// ============================================================================
  
  IXFFDihedral::IXFFDihedral(DihedralType type, int_ id,
                                     std::initializer_list<float_> parameters)
  : _type(type), _id(id) {
    // Allow phase, force constant, multiplicty
    if (_type == DihedralType::Proper) _mask.from_uint32(7);
    // Allow force constant, ideal angle
    if (_type == DihedralType::Improper) _mask.from_uint32(10);
    
    size_ expected_size = _mask.count();
    if (parameters.size() != expected_size)
      throw std::range_error("Incorrect item count");
    
    _dat.fill(0);
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  test_case("IXFFDihedral") {
    using FTyp = DihedralType;
    using TDhd = test::TestFFDihedral;
    using AllowEnum = indigox::test::TestFFDihedral::AllowEnum;
    using StoreEnum = indigox::test::TestFFDihedral::StoreEnum;
    
    subcase("Empty check") {
      check_nothrow(TDhd t(FTyp::Empty, 4, {}));
      check_throws_as(TDhd t(FTyp::Empty, 4, {1.,2.,}), std::range_error);
      
      TDhd testdhd(FTyp::Empty, 4, {});
      check(testdhd.get_mask().none());
      for (size_ i = 0; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testdhd.get_dat()[i]);
      check_throws_as(testdhd.GetPhaseShift(), std::runtime_error);
      check_throws_as(testdhd.GetForceConstant(), std::runtime_error);
      check_throws_as(testdhd.GetIdealAngle(), std::runtime_error);
      check_throws_as(testdhd.GetMultiplicity(), std::runtime_error);
    }
    
    subcase("Proper dihedral check") {
      // throwing construction checks
      check_nothrow(TDhd t(FTyp::Proper, 12, {1.,2.,3.}));
      check_throws_as(TDhd t(FTyp::Proper, 12, {1.,}), std::range_error);
      check_throws_as(TDhd t(FTyp::Proper, 12, {1.,2.,3.,4.}), std::range_error);
      
      TDhd testdhd(FTyp::Proper, 12, {0.00005, 23.09, 5.00007});
      // Mask checks
      check(testdhd.get_mask()[AllowEnum::Allow_PhaseShift]);
      check(testdhd.get_mask()[AllowEnum::Allow_ForceConstant]);
      check(testdhd.get_mask()[AllowEnum::Allow_Multiplicity]);
      check_false(testdhd.get_mask()[AllowEnum::Allow_IdealAngle]);
      // Storage checks
      check_eq(approximately(0.00005), testdhd.get_dat()[0]);
      check_eq(approximately(23.09), testdhd.get_dat()[1]);
      check_eq(approximately(5.00007), testdhd.get_dat()[2]);
      for (size_ i = 3; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testdhd.get_dat()[i]);
      check_eq(12, testdhd.get_id());
      check_eq(FTyp::Proper, testdhd.get_type());
      // Throwing checks
      check_nothrow(testdhd.GetPhaseShift());
      check_nothrow(testdhd.GetForceConstant());
      check_nothrow(testdhd.GetMultiplicity());
      check_throws_as(testdhd.GetIdealAngle(), std::runtime_error);
      // Value checks
      check_eq(approximately(0.00005), testdhd.GetPhaseShift());
      check_eq(approximately(23.09), testdhd.GetForceConstant());
      check_eq(5, testdhd.GetMultiplicity());
      check_eq(12, testdhd.GetID());
      check_eq(FTyp::Proper, testdhd.GetType());
    }
    
    subcase("Improper dihedral check") {
      // throwing construction checks
      check_nothrow(TDhd t(FTyp::Improper, 12, {1.,2.}));
      check_throws_as(TDhd t(FTyp::Improper, 12, {1.,}), std::range_error);
      check_throws_as(TDhd t(FTyp::Improper, 12, {1.,2.,3.}), std::range_error);
      
      TDhd testdhd(FTyp::Improper, 12, {0.00005, 23.09});
      // Mask checks
      check_false(testdhd.get_mask()[AllowEnum::Allow_PhaseShift]);
      check(testdhd.get_mask()[AllowEnum::Allow_ForceConstant]);
      check_false(testdhd.get_mask()[AllowEnum::Allow_Multiplicity]);
      check(testdhd.get_mask()[AllowEnum::Allow_IdealAngle]);
      // Storage checks
      check_eq(approximately(0.00005), testdhd.get_dat()[0]);
      check_eq(approximately(23.09), testdhd.get_dat()[1]);
      for (size_ i = 2; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testdhd.get_dat()[i]);
      check_eq(12, testdhd.get_id());
      check_eq(FTyp::Improper, testdhd.get_type());
      // Throwing checks
      check_throws_as(testdhd.GetPhaseShift(), std::runtime_error);
      check_nothrow(testdhd.GetForceConstant());
      check_throws_as(testdhd.GetMultiplicity(), std::runtime_error);
      check_nothrow(testdhd.GetIdealAngle());
      // Value checks
      check_eq(approximately(0.00005), testdhd.GetIdealAngle());
      check_eq(approximately(23.09), testdhd.GetForceConstant());
      check_eq(12, testdhd.GetID());
      check_eq(FTyp::Improper, testdhd.GetType());
    }
  }
  
// ============================================================================
// == IXFFAngle Serialisation =============================================
// ============================================================================
  
  template <class Archive>
  void IXFFAngle::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("mask", _mask));
  }
  
  template <class Archive>
  void IXFFAngle::load_and_construct(Archive &archive,
                                            cereal::construct<IXFFAngle> &construct,
                                            const uint32_t) {
    construct(AngleType::Empty, 0, std::initializer_list<float_>());
    archive(INDIGOX_SERIAL_NVP("type", construct->_type),
            INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("mask", construct->_mask));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFAngle);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAngle serilisation", T, ixffang_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestFFAngle saved(AngleType::Harmonic, 2, {0.1,0.2});
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestFFAngle loaded(AngleType::Empty, 0, {});
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded.imp));
    }
    
    check_eq(saved.get_id(), loaded.get_id());
    check_eq(saved.get_dat(), loaded.get_dat());
    check_eq(saved.get_mask(), loaded.get_mask());
    check_eq(saved.get_type(), loaded.get_type());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffang_serial, ixserial<IXFFAngle>);
  
// ============================================================================
// == IXFFAngle Construction ==============================================
// ============================================================================
  
  IXFFAngle::IXFFAngle(AngleType type, int_ id,
                               std::initializer_list<float_> parameters)
  : _type(type), _id(id) {
    // Allow force constant and ideal angle
    if (_type == AngleType::Harmonic) _mask.from_uint32(3);
    // Allow force constant and ideal angle
    if (_type == AngleType::CosineHarmonic) _mask.from_uint32(3);
    
    size_ expected_size = _mask.count();
    if (parameters.size() != expected_size)
      throw std::range_error("Incorrect item count");
    
    _dat.fill(0);
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  test_case("IXFFAngle") {
    using FTyp = AngleType;
    using TAng = test::TestFFAngle;
    using AllowEnum = indigox::test::TestFFAngle::AllowEnum;
    using StoreEnum = indigox::test::TestFFAngle::StoreEnum;
    
    subcase("Empty check") {
      check_nothrow(TAng t(FTyp::Empty, 4, {}));
      check_throws_as(TAng t(FTyp::Empty, 4, {1.,2.,}), std::range_error);
      
      TAng testang(FTyp::Empty, 4, {});
      check(testang.get_mask().none());
      for (size_ i = 0; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testang.get_dat()[i]);
      check_throws_as(testang.GetIdealAngle(), std::runtime_error);
      check_throws_as(testang.GetForceConstant(), std::runtime_error);
    }
    
    subcase("Harmonic angle check") {
      // throwing construction checks
      check_nothrow(TAng t(FTyp::Harmonic, 11, {1.,2.}));
      check_throws_as(TAng t(FTyp::Harmonic, 11, {1.,}), std::range_error);
      check_throws_as(TAng t(FTyp::Harmonic, 11, {1.,2.,3.}), std::range_error);
      
      TAng testang(FTyp::Harmonic, 3, {1.3e6, 120.00005});
      // Mask checks
      check(testang.get_mask()[AllowEnum::Allow_ForceConstant]);
      check(testang.get_mask()[AllowEnum::Allow_IdealAngle]);
      // Storage checks
      check_eq(approximately(1.3e6), testang.get_dat()[0]);
      check_eq(approximately(120.00005), testang.get_dat()[1]);
      for (size_ i = 2; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testang.get_dat()[i]);
      check_eq(3, testang.get_id());
      check_eq(FTyp::Harmonic, testang.get_type());
      // Throwing checks
      check_nothrow(testang.GetForceConstant());
      check_nothrow(testang.GetIdealAngle());
      // Value checks
      check_eq(approximately(1.3e6), testang.GetForceConstant());
      check_eq(approximately(120.00005), testang.GetIdealAngle());
      check_eq(3, testang.GetID());
      check_eq(FTyp::Harmonic, testang.GetType());
    }
    
    subcase("Cosine-harmonic angle check") {
      // throwing construction checks
      check_nothrow(TAng t(FTyp::CosineHarmonic, 11, {1.,2.}));
      check_throws_as(TAng t(FTyp::CosineHarmonic, 11, {1.,}), std::range_error);
      check_throws_as(TAng t(FTyp::CosineHarmonic, 11, {1.,2.,3.}), std::range_error);
      
      TAng testang(FTyp::CosineHarmonic, 3, {1.3e6, 120.00005});
      // Mask checks
      check(testang.get_mask()[AllowEnum::Allow_ForceConstant]);
      check(testang.get_mask()[AllowEnum::Allow_IdealAngle]);
      // Storage checks
      check_eq(approximately(1.3e6), testang.get_dat()[0]);
      check_eq(approximately(120.00005), testang.get_dat()[1]);
      for (size_ i = 2; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testang.get_dat()[i]);
      check_eq(3, testang.get_id());
      check_eq(FTyp::CosineHarmonic, testang.get_type());
      // Throwing checks
      check_nothrow(testang.GetForceConstant());
      check_nothrow(testang.GetIdealAngle());
      // Value checks
      check_eq(approximately(1.3e6), testang.GetForceConstant());
      check_eq(approximately(120.00005), testang.GetIdealAngle());
      check_eq(3, testang.GetID());
      check_eq(FTyp::CosineHarmonic, testang.GetType());
    }
  }

// ============================================================================
// == IXFFBond Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXFFBond::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("mask", _mask));
  }
  
  template <class Archive>
  void IXFFBond::load_and_construct(Archive &archive,
                                        cereal::construct<IXFFBond> &construct,
                                        const uint32_t) {
    construct(BondType::Empty, 0, std::initializer_list<float_>());
    archive(INDIGOX_SERIAL_NVP("type", construct->_type),
            INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("mask", construct->_mask));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFBond);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFBond serilisation", T, ixffbnd_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestFFBond saved(BondType::Harmonic, 2, {0.1,0.2});
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestFFBond loaded(BondType::Empty, 0, {});
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded.imp));
    }
    
    check_eq(saved.get_id(), loaded.get_id());
    check_eq(saved.get_dat(), loaded.get_dat());
    check_eq(saved.get_mask(), loaded.get_mask());
    check_eq(saved.get_type(), loaded.get_type());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffbnd_serial, ixserial<IXFFBond>);
  
// ============================================================================
// == IXFFBond Construction ===============================================
// ============================================================================
  
  IXFFBond::IXFFBond(BondType type, int_ id,
                             std::initializer_list<float_> parameters)
  : _type(type), _id(id) {
    // Allow force constant and ideal length
    if (_type == BondType::Harmonic) _mask.from_uint32(3);
    // Allow force constant and ideal length
    if (_type == BondType::Quartic) _mask.from_uint32(3);
    
    size_ expected_size = _mask.count();
    if (parameters.size() != expected_size)
      throw std::range_error("Incorrect item count");
    
    _dat.fill(0);
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  test_case("IXFFBond") {
    using FTyp = BondType;
    using TBnd = test::TestFFBond;
    using AllowEnum = indigox::test::TestFFBond::AllowEnum;
    using StoreEnum = indigox::test::TestFFBond::StoreEnum;
    
    subcase("Empty check") {
      check_nothrow(TBnd t(FTyp::Empty, 4, {}));
      check_throws_as(TBnd t(FTyp::Empty, 4, {1.,2.,}), std::range_error);
      
      TBnd testbnd(FTyp::Empty, 4, {});
      check(testbnd.get_mask().none());
      for (size_ i = 0; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testbnd.get_dat()[i]);
      check_throws_as(testbnd.GetIdealLength(), std::runtime_error);
      check_throws_as(testbnd.GetForceConstant(), std::runtime_error);
    }
    
    subcase("Harmonic bond check") {
      // throwing construction checks
      check_nothrow(TBnd t(FTyp::Harmonic, 11, {1.,2.}));
      check_throws_as(TBnd t(FTyp::Harmonic, 11, {1.,}), std::range_error);
      check_throws_as(TBnd t(FTyp::Harmonic, 11, {1.,2.,3.}), std::range_error);
      
      TBnd testbnd(FTyp::Harmonic, 3, {1.3e6, 120.00005});
      // Mask checks
      check(testbnd.get_mask()[AllowEnum::Allow_ForceConstant]);
      check(testbnd.get_mask()[AllowEnum::Allow_IdealLength]);
      // Storage checks
      check_eq(approximately(1.3e6), testbnd.get_dat()[0]);
      check_eq(approximately(120.00005), testbnd.get_dat()[1]);
      for (size_ i = 2; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testbnd.get_dat()[i]);
      check_eq(3, testbnd.get_id());
      check_eq(FTyp::Harmonic, testbnd.get_type());
      // Throwing checks
      check_nothrow(testbnd.GetForceConstant());
      check_nothrow(testbnd.GetIdealLength());
      // Value checks
      check_eq(approximately(1.3e6), testbnd.GetForceConstant());
      check_eq(approximately(120.00005), testbnd.GetIdealLength());
      check_eq(3, testbnd.GetID());
      check_eq(FTyp::Harmonic, testbnd.GetType());
    }
    
    subcase("Quartic bond check") {
      // throwing construction checks
      check_nothrow(TBnd t(FTyp::Quartic, 11, {1.,2.}));
      check_throws_as(TBnd t(FTyp::Quartic, 11, {1.,}), std::range_error);
      check_throws_as(TBnd t(FTyp::Quartic, 11, {1.,2.,3.}), std::range_error);
      
      TBnd testbnd(FTyp::Quartic, 3, {1.3e6, 120.00005});
      // Mask checks
      check(testbnd.get_mask()[AllowEnum::Allow_ForceConstant]);
      check(testbnd.get_mask()[AllowEnum::Allow_IdealLength]);
      // Storage checks
      check_eq(approximately(1.3e6), testbnd.get_dat()[0]);
      check_eq(approximately(120.00005), testbnd.get_dat()[1]);
      for (size_ i = 2; i < StoreEnum::STORE_SIZE; ++i)
        check_eq(approximately(0.0), testbnd.get_dat()[i]);
      check_eq(3, testbnd.get_id());
      check_eq(FTyp::Quartic, testbnd.get_type());
      // Throwing checks
      check_nothrow(testbnd.GetForceConstant());
      check_nothrow(testbnd.GetIdealLength());
      // Value checks
      check_eq(approximately(1.3e6), testbnd.GetForceConstant());
      check_eq(approximately(120.00005), testbnd.GetIdealLength());
      check_eq(3, testbnd.GetID());
      check_eq(FTyp::Quartic, testbnd.GetType());
    }
  }
  
// ============================================================================
// == IXFFAtom Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXFFAtom::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("name", _name));
  }
  
  template <class Archive>
  void IXFFAtom::load_and_construct(Archive &archive,
                                        cereal::construct<IXFFAtom> &construct,
                                        const uint32_t) {
    construct(0, "blank");
    archive(INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("name", construct->_name));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFAtom);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAtom serialisation", T, ixffatm_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestFFAtom saved(12, "tester");
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestFFAtom loaded(0,"");
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded.imp));
    }
    
    check_eq(saved.get_id(), loaded.get_id());
    check_eq(saved.get_name(), loaded.get_name());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffatm_serial, ixserial<IXFFAtom>);
  
// ============================================================================
// == IXFFAtom Construction ===============================================
// ============================================================================
  
  IXFFAtom::IXFFAtom(int_ id, string_ name) : _id(id), _name(name) { }
  
  test_case("IXFFAtom") {
    subcase("Construction check") {
      check_nothrow(test::TestFFAtom t(0, ""));
      
      test::TestFFAtom testatm(2, "CH0");
      check_eq(2, testatm.GetID());
      check_eq("CH0", testatm.GetName());
    }
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
  
  IXForcefield::IXForcefield(FFFamily family, string_ name)
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
  
  FFAtom IXForcefield::NewAtomType(int_ id, string_ name) {
    if (GetAtomType(id) || GetAtomType(name)) return FFAtom();
    _atms.emplace_back(new IXFFAtom(id, name));
    return _atms.back();
  }
  
  FFBond IXForcefield::NewBondType(BondType type, int_ id,
                                       std::initializer_list<float_> parameters) {
    if (_bnds.find(type) == _bnds.end()) return FFBond();
    if (GetBondType(type, id)) return FFBond();
    _bnds[type].emplace_back(new IXFFBond(type, id, parameters));
    return _bnds[type].back();
  }
  
  FFAngle IXForcefield::NewAngleType(AngleType type, int_ id,
                                       std::initializer_list<float_> parameters) {
    if (_angs.find(type) == _angs.end()) return FFAngle();
    if (GetAngleType(type, id)) return FFAngle();
    _angs[type].emplace_back(new IXFFAngle(type, id, parameters));
    return _angs[type].back();
  }
  
  FFDihedral IXForcefield::NewDihedralType(DihedralType type, int_ id,
                                       std::initializer_list<float_> parameters) {
    if (_dhds.find(type) == _dhds.end()) return FFDihedral();
    if (GetDihedralType(type, id)) return FFDihedral();
    _dhds[type].emplace_back(new IXFFDihedral(type, id, parameters));
    return _dhds[type].back();
  }
  
// ============================================================================
// == IXForcefield Getting types ==============================================
// ============================================================================
  
  FFAtom IXForcefield::GetAtomType(string_ name) const {
    auto pred = [&name](const FFAtom& t) { return t->GetName() == name; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? FFAtom() : *pos;
  }
  
  FFAtom IXForcefield::GetAtomType(int_ id) const {
    auto pred = [&id](const FFAtom& t) { return t->GetID() == id; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? FFAtom() : *pos;
  }
  
  FFBond IXForcefield::GetBondType(BondType type, int_ id) const {
    if (_bnds.find(type) == _bnds.end()) return FFBond();
    auto pred = [&id](const FFBond& t) { return t->GetID() == id; };
    auto pos = std::find_if(_bnds.at(type).begin(), _bnds.at(type).end(), pred);
    return (pos == _bnds.at(type).end()) ? FFBond() : *pos;
  }
  
  FFAngle IXForcefield::GetAngleType(AngleType type, int_ id) const {
    if (_angs.find(type) == _angs.end()) return FFAngle();
    auto pred = [&id](const FFAngle& t) { return t->GetID() == id; };
    auto pos = std::find_if(_angs.at(type).begin(), _angs.at(type).end(), pred);
    return (pos == _angs.at(type).end()) ? FFAngle() : *pos;
  }
  
  FFDihedral IXForcefield::GetDihedralType(DihedralType type, int_ id) const {
    if (_dhds.find(type) == _dhds.end()) return FFDihedral();
    auto pred = [&id](const FFDihedral& t) { return t->GetID() == id; };
    auto pos = std::find_if(_dhds.at(type).begin(), _dhds.at(type).end(), pred);
    return (pos == _dhds.at(type).end()) ? FFDihedral() : *pos;
  }
  
  test_case("IXForcefield") {
    test::TestForcefield ff(FFFamily::GROMOS, "Test");
    ff.imp = GenerateGROMOS54A7();
    check_eq(54, ff.NumAtomTypes());
    check_eq(104, ff.NumBondTypes());
    check_eq(108, ff.NumAngleTypes());
    check_eq(50, ff.NumDihedralTypes());
  }
  
  Forcefield GenerateGROMOS54A7() {
    Forcefield ff = std::make_shared<IXForcefield>(FFFamily::GROMOS, "54A7");
    // Add bond types
    ff->NewHarmonicBondType(1, 3.1400000e+05, 1.0000000e-01);
    ff->NewQuarticBondType(1, 1.5700000e+07, 1.0000000e-01);
    ff->NewHarmonicBondType(2, 3.7400000e+05, 1.0000000e-01);
    ff->NewQuarticBondType(2, 1.8700000e+07, 1.0000000e-01);
    ff->NewHarmonicBondType(3, 2.9227260e+05, 1.0900000e-01);
    ff->NewQuarticBondType(3, 1.2300000e+07, 1.0900000e-01);
    ff->NewHarmonicBondType(4, 9.2825600e+05, 1.1200000e-01);
    ff->NewQuarticBondType(4, 3.7000000e+07, 1.1200000e-01);
    ff->NewHarmonicBondType(5, 5.0228280e+05, 1.2300000e-01);
    ff->NewQuarticBondType(5, 1.6600000e+07, 1.2300000e-01);
    ff->NewHarmonicBondType(6, 4.1875000e+05, 1.2500000e-01);
    ff->NewQuarticBondType(6, 1.3400000e+07, 1.2500000e-01);
    ff->NewHarmonicBondType(7, 4.1817600e+05, 1.3200000e-01);
    ff->NewQuarticBondType(7, 1.2000000e+07, 1.3200000e-01);
    ff->NewHarmonicBondType(8, 3.1380286e+05, 1.3300000e-01);
    ff->NewQuarticBondType(8, 8.8700000e+06, 1.3300000e-01);
    ff->NewHarmonicBondType(9, 3.7500680e+05, 1.3300000e-01);
    ff->NewQuarticBondType(9, 1.0600000e+07, 1.3300000e-01);
    ff->NewHarmonicBondType(10, 4.1746040e+05, 1.3300000e-01);
    ff->NewQuarticBondType(10, 1.1800000e+07, 1.3300000e-01);
    ff->NewHarmonicBondType(11, 3.7707600e+05, 1.3400000e-01);
    ff->NewQuarticBondType(11, 1.0500000e+07, 1.3400000e-01);
    ff->NewHarmonicBondType(12, 4.2017040e+05, 1.3400000e-01);
    ff->NewQuarticBondType(12, 1.1700000e+07, 1.3400000e-01);
    ff->NewHarmonicBondType(13, 3.7731840e+05, 1.3600000e-01);
    ff->NewQuarticBondType(13, 1.0200000e+07, 1.3600000e-01);
    ff->NewHarmonicBondType(14, 4.1896800e+05, 1.3800000e-01);
    ff->NewQuarticBondType(14, 1.1000000e+07, 1.3800000e-01);
    ff->NewHarmonicBondType(15, 3.3463972e+05, 1.3900000e-01);
    ff->NewQuarticBondType(15, 8.6600000e+06, 1.3900000e-01);
    ff->NewHarmonicBondType(16, 4.1733360e+05, 1.3900000e-01);
    ff->NewQuarticBondType(16, 1.0800000e+07, 1.3900000e-01);
    ff->NewHarmonicBondType(17, 3.3476800e+05, 1.4000000e-01);
    ff->NewQuarticBondType(17, 8.5400000e+06, 1.4000000e-01);
    ff->NewHarmonicBondType(18, 3.3454564e+05, 1.4300000e-01);
    ff->NewQuarticBondType(18, 8.1800000e+06, 1.4300000e-01);
    ff->NewHarmonicBondType(19, 3.7667058e+05, 1.4300000e-01);
    ff->NewQuarticBondType(19, 9.2100000e+06, 1.4300000e-01);
    ff->NewHarmonicBondType(20, 2.5122545e+05, 1.4350000e-01);
    ff->NewQuarticBondType(20, 6.1000000e+06, 1.4350000e-01);
    ff->NewHarmonicBondType(21, 3.7642878e+05, 1.4700000e-01);
    ff->NewQuarticBondType(21, 8.7100000e+06, 1.4700000e-01);
    ff->NewHarmonicBondType(22, 2.5101984e+05, 1.4800000e-01);
    ff->NewQuarticBondType(22, 5.7300000e+06, 1.4800000e-01);
    ff->NewHarmonicBondType(23, 3.3469312e+05, 1.4800000e-01);
    ff->NewQuarticBondType(23, 7.6400000e+06, 1.4800000e-01);
    ff->NewHarmonicBondType(24, 3.7674880e+05, 1.4800000e-01);
    ff->NewQuarticBondType(24, 8.6000000e+06, 1.4800000e-01);
    ff->NewHarmonicBondType(25, 3.7665000e+05, 1.5000000e-01);
    ff->NewQuarticBondType(25, 8.3700000e+06, 1.5000000e-01);
    ff->NewHarmonicBondType(26, 2.5090944e+05, 1.5200000e-01);
    ff->NewQuarticBondType(26, 5.4300000e+06, 1.5200000e-01);
    ff->NewHarmonicBondType(27, 3.3474870e+05, 1.5300000e-01);
    ff->NewQuarticBondType(27, 7.1500000e+06, 1.5300000e-01);
    ff->NewHarmonicBondType(28, 2.5091528e+05, 1.6100000e-01);
    ff->NewQuarticBondType(28, 4.8400000e+06, 1.6100000e-01);
    ff->NewHarmonicBondType(29, 2.5081136e+05, 1.6300000e-01);
    ff->NewQuarticBondType(29, 4.7200000e+06, 1.6300000e-01);
    ff->NewHarmonicBondType(30, 1.7236096e+05, 1.7800000e-01);
    ff->NewQuarticBondType(30, 2.7200000e+06, 1.7800000e-01);
    ff->NewHarmonicBondType(31, 3.7640592e+05, 1.7800000e-01);
    ff->NewQuarticBondType(31, 5.9400000e+06, 1.7800000e-01);
    ff->NewHarmonicBondType(32, 3.7641636e+05, 1.8300000e-01);
    ff->NewQuarticBondType(32, 5.6200000e+06, 1.8300000e-01);
    ff->NewHarmonicBondType(33, 2.5107742e+05, 1.8700000e-01);
    ff->NewQuarticBondType(33, 3.5900000e+06, 1.8700000e-01);
    ff->NewHarmonicBondType(34, 5.0181120e+04, 1.9800000e-01);
    ff->NewQuarticBondType(34, 6.4000000e+05, 1.9800000e-01);
    ff->NewHarmonicBondType(35, 5.0240000e+04, 2.0000000e-01);
    ff->NewQuarticBondType(35, 6.2800000e+05, 2.0000000e-01);
    ff->NewHarmonicBondType(36, 4.1865696e+05, 2.0400000e-01);
    ff->NewQuarticBondType(36, 5.0300000e+06, 2.0400000e-01);
    ff->NewHarmonicBondType(37, 5.2748280e+04, 2.2100000e-01);
    ff->NewQuarticBondType(37, 5.4000000e+05, 2.2100000e-01);
    ff->NewHarmonicBondType(38, 4.6400000e+05, 1.0000000e-01);
    ff->NewQuarticBondType(38, 2.3200000e+07, 1.0000000e-01);
    ff->NewHarmonicBondType(39, 2.9282000e+05, 1.1000000e-01);
    ff->NewQuarticBondType(39, 1.2100000e+07, 1.1000000e-01);
    ff->NewHarmonicBondType(40, 5.0190759e+05, 1.7580000e-01);
    ff->NewQuarticBondType(40, 8.1200000e+06, 1.7580000e-01);
    ff->NewHarmonicBondType(41, 3.7641672e+05, 1.5300000e-01);
    ff->NewQuarticBondType(41, 8.0400000e+06, 1.5300000e-01);
    ff->NewHarmonicBondType(42, 3.7182472e+05, 1.9379900e-01);
    ff->NewQuarticBondType(42, 4.9500000e+06, 1.9379900e-01);
    ff->NewHarmonicBondType(43, 5.0181120e+05, 1.7600000e-01);
    ff->NewQuarticBondType(43, 8.1000000e+06, 1.7600000e-01);
    ff->NewHarmonicBondType(44, 4.1925895e+05, 1.2650000e-01);
    ff->NewQuarticBondType(44, 1.3100000e+07, 1.2650000e-01);
    ff->NewHarmonicBondType(45, 3.7543500e+05, 1.3500000e-01);
    ff->NewQuarticBondType(45, 1.0300000e+07, 1.3500000e-01);
    ff->NewHarmonicBondType(46, 4.6453153e+05, 1.6329900e-01);
    ff->NewQuarticBondType(46, 8.7100000e+06, 1.6329900e-01);
    ff->NewHarmonicBondType(47, 2.9308843e+05, 2.3383900e-01);
    ff->NewQuarticBondType(47, 2.6800000e+06, 2.3383900e-01);
    ff->NewHarmonicBondType(48, 5.0221475e+05, 2.9028300e-01);
    ff->NewQuarticBondType(48, 2.9800000e+06, 2.9028300e-01);
    ff->NewHarmonicBondType(49, 3.7311559e+05, 2.7938800e-01);
    ff->NewQuarticBondType(49, 2.3900000e+06, 2.7938800e-01);
    ff->NewHarmonicBondType(50, 3.7138473e+05, 2.9118900e-01);
    ff->NewQuarticBondType(50, 2.1900000e+06, 2.9118900e-01);
    ff->NewHarmonicBondType(51, 3.4252596e+05, 2.0770000e-01);
    ff->NewQuarticBondType(51, 3.9700000e+06, 2.0770000e-01);
    ff->NewHarmonicBondType(52, 5.0222492e+05, 2.8740700e-01);
    ff->NewQuarticBondType(52, 3.0400000e+06, 2.8740700e-01);
    for (int_ i = 1; i < 53; ++i)
      ff->LinkBondTypes(ff->GetHarmonicBondType(i), ff->GetQuarticBondType(i));
    
    
    // Add Angle types
    ff->NewHarmonicAngleType(1, 1.1550101e-01, 9.0000000e+01);
    ff->NewCosineHarmonicAngleType(1, 3.8000000e+02, 9.0000000e+01);
    ff->NewHarmonicAngleType(2, 1.2768574e-01, 9.0000000e+01);
    ff->NewCosineHarmonicAngleType(2, 4.2000000e+02, 9.0000000e+01);
    ff->NewHarmonicAngleType(3, 1.2177061e-01, 9.6000000e+01);
    ff->NewCosineHarmonicAngleType(3, 4.0500000e+02, 9.6000000e+01);
    ff->NewHarmonicAngleType(4, 1.4008261e-01, 1.0000000e+02);
    ff->NewCosineHarmonicAngleType(4, 4.7500000e+02, 1.0000000e+02);
    ff->NewHarmonicAngleType(5, 1.2122177e-01, 1.0300000e+02);
    ff->NewCosineHarmonicAngleType(5, 4.2000000e+02, 1.0300000e+02);
    ff->NewHarmonicAngleType(6, 1.4028506e-01, 1.0400000e+02);
    ff->NewCosineHarmonicAngleType(6, 4.9000000e+02, 1.0400000e+02);
    ff->NewHarmonicAngleType(7, 1.2788754e-01, 1.0800000e+02);
    ff->NewCosineHarmonicAngleType(7, 4.6500000e+02, 1.0800000e+02);
    ff->NewHarmonicAngleType(8, 7.6912479e-02, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(8, 2.8500000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(9, 8.6386140e-02, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(9, 3.2000000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(10, 1.0262668e-01, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(10, 3.8000000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(11, 1.1480708e-01, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(11, 4.2500000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(12, 1.2157397e-01, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(12, 4.5000000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(13, 1.4052124e-01, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(13, 5.2000000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(14, 1.2142334e-01, 1.0960000e+02);
    ff->NewCosineHarmonicAngleType(14, 4.5000000e+02, 1.0960000e+02);
    ff->NewHarmonicAngleType(15, 1.4048747e-01, 1.1100000e+02);
    ff->NewCosineHarmonicAngleType(15, 5.3000000e+02, 1.1100000e+02);
    ff->NewHarmonicAngleType(16, 1.4045138e-01, 1.1300000e+02);
    ff->NewCosineHarmonicAngleType(16, 5.4500000e+02, 1.1300000e+02);
    ff->NewHarmonicAngleType(17, 1.2296570e-02, 1.1500000e+02);
    ff->NewCosineHarmonicAngleType(17, 5.0000000e+01, 1.1500000e+02);
    ff->NewHarmonicAngleType(18, 1.1488482e-01, 1.1500000e+02);
    ff->NewCosineHarmonicAngleType(18, 4.6000000e+02, 1.1500000e+02);
    ff->NewHarmonicAngleType(19, 1.5241650e-01, 1.1500000e+02);
    ff->NewCosineHarmonicAngleType(19, 6.1000000e+02, 1.1500000e+02);
    ff->NewHarmonicAngleType(20, 1.1421859e-01, 1.1600000e+02);
    ff->NewCosineHarmonicAngleType(20, 4.6500000e+02, 1.1600000e+02);
    ff->NewHarmonicAngleType(21, 1.5236094e-01, 1.1600000e+02);
    ff->NewCosineHarmonicAngleType(21, 6.2000000e+02, 1.1600000e+02);
    ff->NewHarmonicAngleType(22, 1.5336019e-01, 1.1700000e+02);
    ff->NewCosineHarmonicAngleType(22, 6.3500000e+02, 1.1700000e+02);
    ff->NewHarmonicAngleType(23, 8.8910434e-02, 1.2000000e+02);
    ff->NewCosineHarmonicAngleType(23, 3.9000000e+02, 1.2000000e+02);
    ff->NewHarmonicAngleType(24, 1.0147593e-01, 1.2000000e+02);
    ff->NewCosineHarmonicAngleType(24, 4.4500000e+02, 1.2000000e+02);
    ff->NewHarmonicAngleType(25, 1.1518373e-01, 1.2000000e+02);
    ff->NewCosineHarmonicAngleType(25, 5.0500000e+02, 1.2000000e+02);
    ff->NewHarmonicAngleType(26, 1.2089532e-01, 1.2000000e+02);
    ff->NewCosineHarmonicAngleType(26, 5.3000000e+02, 1.2000000e+02);
    ff->NewHarmonicAngleType(27, 1.2774922e-01, 1.2000000e+02);
    ff->NewCosineHarmonicAngleType(27, 5.6000000e+02, 1.2000000e+02);
    ff->NewHarmonicAngleType(28, 1.5288018e-01, 1.2000000e+02);
    ff->NewCosineHarmonicAngleType(28, 6.7000000e+02, 1.2000000e+02);
    ff->NewHarmonicAngleType(29, 1.7801113e-01, 1.2000000e+02);
    ff->NewCosineHarmonicAngleType(29, 7.8000000e+02, 1.2000000e+02);
    ff->NewHarmonicAngleType(30, 1.5312732e-01, 1.2100000e+02);
    ff->NewCosineHarmonicAngleType(30, 6.8500000e+02, 1.2100000e+02);
    ff->NewHarmonicAngleType(31, 1.5317431e-01, 1.2200000e+02);
    ff->NewCosineHarmonicAngleType(31, 7.0000000e+02, 1.2200000e+02);
    ff->NewHarmonicAngleType(32, 8.8743846e-02, 1.2300000e+02);
    ff->NewCosineHarmonicAngleType(32, 4.1500000e+02, 1.2300000e+02);
    ff->NewHarmonicAngleType(33, 1.5266919e-01, 1.2400000e+02);
    ff->NewCosineHarmonicAngleType(33, 7.3000000e+02, 1.2400000e+02);
    ff->NewHarmonicAngleType(34, 7.6490216e-02, 1.2500000e+02);
    ff->NewCosineHarmonicAngleType(34, 3.7500000e+02, 1.2500000e+02);
    ff->NewHarmonicAngleType(35, 1.5314080e-01, 1.2500000e+02);
    ff->NewCosineHarmonicAngleType(35, 7.5000000e+02, 1.2500000e+02);
    ff->NewHarmonicAngleType(36, 1.1448735e-01, 1.2600000e+02);
    ff->NewCosineHarmonicAngleType(36, 5.7500000e+02, 1.2600000e+02);
    ff->NewHarmonicAngleType(37, 1.2744672e-01, 1.2600000e+02);
    ff->NewCosineHarmonicAngleType(37, 6.4000000e+02, 1.2600000e+02);
    ff->NewHarmonicAngleType(38, 1.5336544e-01, 1.2600000e+02);
    ff->NewCosineHarmonicAngleType(38, 7.7000000e+02, 1.2600000e+02);
    ff->NewHarmonicAngleType(39, 1.2775497e-01, 1.3200000e+02);
    ff->NewCosineHarmonicAngleType(39, 7.6000000e+02, 1.3200000e+02);
    ff->NewHarmonicAngleType(40, 1.2112698e-01, 1.5500000e+02);
    ff->NewCosineHarmonicAngleType(40, 2.2150000e+03, 1.5500000e+02);
    ff->NewHarmonicAngleType(41, 7.2640156e-02, 1.8000000e+02);
    ff->NewCosineHarmonicAngleType(41, 9.1350000e+04, 1.8000000e+02);
    ff->NewHarmonicAngleType(42, 1.1724316e-01, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(42, 4.3400000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(43, 1.3376523e-01, 1.0757000e+02);
    ff->NewCosineHarmonicAngleType(43, 4.8400000e+02, 1.0757000e+02);
    ff->NewHarmonicAngleType(44, 1.6689058e-01, 1.1130000e+02);
    ff->NewCosineHarmonicAngleType(44, 6.3200000e+02, 1.1130000e+02);
    ff->NewHarmonicAngleType(45, 1.4024534e-01, 9.7400000e+01);
    ff->NewCosineHarmonicAngleType(45, 4.6900000e+02, 9.7400000e+01);
    ff->NewHarmonicAngleType(46, 1.4026005e-01, 1.0675000e+02);
    ff->NewCosineHarmonicAngleType(46, 5.0300000e+02, 1.0675000e+02);
    ff->NewHarmonicAngleType(47, 1.2108416e-01, 1.0853000e+02);
    ff->NewCosineHarmonicAngleType(47, 4.4300000e+02, 1.0853000e+02);
    ff->NewHarmonicAngleType(48, 1.6704740e-01, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(48, 6.1800000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(49, 1.4008648e-01, 1.0760000e+02);
    ff->NewCosineHarmonicAngleType(49, 5.0700000e+02, 1.0760000e+02);
    ff->NewHarmonicAngleType(50, 1.2103261e-01, 1.0950000e+02);
    ff->NewCosineHarmonicAngleType(50, 4.4800000e+02, 1.0950000e+02);
    ff->NewHarmonicAngleType(51, 1.4017954e-01, 1.1030000e+02);
    ff->NewCosineHarmonicAngleType(51, 5.2400000e+02, 1.1030000e+02);
    ff->NewHarmonicAngleType(52, 1.4025677e-01, 1.1140000e+02);
    ff->NewCosineHarmonicAngleType(52, 5.3200000e+02, 1.1140000e+02);
    ff->NewHarmonicAngleType(53, 1.5305438e-01, 1.1720000e+02);
    ff->NewCosineHarmonicAngleType(53, 6.3600000e+02, 1.1720000e+02);
    ff->NewHarmonicAngleType(54, 1.5294820e-01, 1.2140000e+02);
    ff->NewCosineHarmonicAngleType(54, 6.9000000e+02, 1.2140000e+02);
    for (int_ i = 1; i < 55; ++i)
      ff->LinkAngleTypes(ff->GetHarmonicAngleType(i),
                         ff->GetCosineHarmonicAngleType(i));
    
    // Add Improper types
    ff->NewImproperDihedralType(1, 0.0510, 0.0);
    ff->NewImproperDihedralType(2, 0.102, 35.26439);
    ff->NewImproperDihedralType(3, 0.204, 0.0);
    ff->NewImproperDihedralType(4, 0.0510, 180.0);
    ff->NewImproperDihedralType(5, 0.102, -35.26439);
    
    // Add Proper types
    ff->NewProperDihedralType(1, 2.670, 180.0, 1);
    ff->NewProperDihedralType(2, 3.410, 180.0, 1);
    ff->NewProperDihedralType(3, 4.970, 180.0, 1);
    ff->NewProperDihedralType(4, 5.860, 180.0, 1);
    ff->NewProperDihedralType(5, 9.350, 180.0, 1);
    ff->NewProperDihedralType(6, 9.450, 180.0, 1);
    ff->NewProperDihedralType(7, 2.790, 0.0, 1);
    ff->NewProperDihedralType(8, 5.350, 0.0, 1);
    ff->NewProperDihedralType(9, 1.530, 180.0, 2);
    ff->NewProperDihedralType(10, 5.860, 180.0, 2);
    ff->NewProperDihedralType(11, 7.110, 180.0, 2);
    ff->NewProperDihedralType(12, 16.700, 180.0, 2);
    ff->NewProperDihedralType(13, 24.000, 180.0, 2);
    ff->NewProperDihedralType(14, 33.500, 180.0, 2);
    ff->NewProperDihedralType(15, 41.800, 180.0, 2);
    ff->NewProperDihedralType(16, 0.000, 0.0, 2);
    ff->NewProperDihedralType(17, 0.418, 0.0, 2);
    ff->NewProperDihedralType(18, 2.090, 0.0, 2);
    ff->NewProperDihedralType(19, 3.140, 0.0, 2);
    ff->NewProperDihedralType(20, 5.090, 0.0, 2);
    ff->NewProperDihedralType(21, 16.700, 0.0, 2);
    ff->NewProperDihedralType(22, 1.050, 0.0, 3);
    ff->NewProperDihedralType(23, 1.260, 0.0, 3);
    ff->NewProperDihedralType(24, 1.300, 0.0, 3);
    ff->NewProperDihedralType(25, 2.530, 0.0, 3);
    ff->NewProperDihedralType(26, 2.930, 0.0, 3);
    ff->NewProperDihedralType(27, 3.190, 0.0, 3);
    ff->NewProperDihedralType(28, 3.650, 0.0, 3);
    ff->NewProperDihedralType(29, 3.770, 0.0, 3);
    ff->NewProperDihedralType(30, 3.900, 0.0, 3);
    ff->NewProperDihedralType(31, 4.180, 0.0, 3);
    ff->NewProperDihedralType(32, 4.690, 0.0, 3);
    ff->NewProperDihedralType(33, 5.440, 0.0, 3);
    ff->NewProperDihedralType(34, 5.920, 0.0, 3);
    ff->NewProperDihedralType(35, 7.690, 0.0, 3);
    ff->NewProperDihedralType(36, 8.620, 0.0, 3);
    ff->NewProperDihedralType(37, 9.500, 0.0, 3);
    ff->NewProperDihedralType(38, 0.000, 0.0, 4);
    ff->NewProperDihedralType(39, 1.000, 180.0, 6);
    ff->NewProperDihedralType(40, 1.000, 0.0, 6);
    ff->NewProperDihedralType(41, 3.770, 0.0, 6);
    ff->NewProperDihedralType(42, 3.500, 180.0, 2);
    ff->NewProperDihedralType(43, 2.800, 0.0, 3);
    ff->NewProperDihedralType(44, 0.700, 180.0, 6);
    ff->NewProperDihedralType(45, 0.400, 0.0, 6);
    
    // Add Atom types
    ff->NewAtomType(1, "O");
    ff->NewAtomType(2, "OM");
    ff->NewAtomType(3, "OA");
    ff->NewAtomType(4, "OE");
    ff->NewAtomType(5, "OW");
    ff->NewAtomType(6, "N");
    ff->NewAtomType(7, "NT");
    ff->NewAtomType(8, "NL");
    ff->NewAtomType(9, "NR");
    ff->NewAtomType(10, "NZ");
    ff->NewAtomType(11, "NE");
    ff->NewAtomType(12, "C");
    ff->NewAtomType(13, "CH0");
    ff->NewAtomType(14, "CH1");
    ff->NewAtomType(15, "CH2");
    ff->NewAtomType(16, "CH3");
    ff->NewAtomType(17, "CH4");
    ff->NewAtomType(18, "CH2r");
    ff->NewAtomType(19, "CR1");
    ff->NewAtomType(20, "HC");
    ff->NewAtomType(21, "H");
    ff->NewAtomType(22, "DUM");
    ff->NewAtomType(23, "S");
    ff->NewAtomType(24, "CU1+");
    ff->NewAtomType(25, "CU2+");
    ff->NewAtomType(26, "FE");
    ff->NewAtomType(27, "ZN2+");
    ff->NewAtomType(28, "MG2+");
    ff->NewAtomType(29, "CA2+");
    ff->NewAtomType(30, "P,SI");
    ff->NewAtomType(31, "AR");
    ff->NewAtomType(32, "F");
    ff->NewAtomType(33, "CL");
    ff->NewAtomType(34, "BR");
    ff->NewAtomType(35, "CMet");
    ff->NewAtomType(36, "OMet");
    ff->NewAtomType(37, "NA+");
    ff->NewAtomType(38, "CL-");
    ff->NewAtomType(39, "CChl");
    ff->NewAtomType(40, "CLChl");
    ff->NewAtomType(41, "HChl");
    ff->NewAtomType(42, "SDmso");
    ff->NewAtomType(43, "CDmso");
    ff->NewAtomType(44, "ODmso");
    ff->NewAtomType(45, "CCl4");
    ff->NewAtomType(46, "CLCl4");
    ff->NewAtomType(47, "FTFE");
    ff->NewAtomType(48, "CTFE");
    ff->NewAtomType(49, "CHTFE");
    ff->NewAtomType(50, "OTFE");
    ff->NewAtomType(51, "CUrea");
    ff->NewAtomType(52, "OUrea");
    ff->NewAtomType(53, "NUrea");
    ff->NewAtomType(54, "CH3p");
    return ff;
  }
  
  test_suite_close();
}
