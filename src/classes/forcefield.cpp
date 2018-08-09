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
// == IXFFDihedralType Serialisation ==========================================
// ============================================================================
  
  template <class Archive>
  void IXFFDihedralType::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("mask", _mask));
  }
  
  template <class Archive>
  void IXFFDihedralType::load_and_construct(Archive &archive,
                                            cereal::construct<IXFFDihedralType> &construct,
                                            const uint32_t) {
    construct(DihedralFunctionType::Empty, 0, std::initializer_list<float_>());
    archive(INDIGOX_SERIAL_NVP("type", construct->_type),
            INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("mask", construct->_mask));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFDihedralType);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFDihedral serilisation", T, ixffdhd_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestDihedralType saved(DihedralFunctionType::Proper, 72, {0.1,0.2,0.3});
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestDihedralType loaded(DihedralFunctionType::Empty, 0, {});
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
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffdhd_serial, ixserial<IXFFDihedralType>);
  
// ============================================================================
// == IXFFDihedralType Construction ===========================================
// ============================================================================
  
  IXFFDihedralType::IXFFDihedralType(DihedralFunctionType type, int_ id,
                                     std::initializer_list<float_> parameters)
  : _type(type), _id(id) {
    // Allow phase, force constant, multiplicty
    if (_type == DihedralFunctionType::Proper) _mask.from_uint32(7);
    // Allow force constant, ideal angle
    if (_type == DihedralFunctionType::Improper) _mask.from_uint32(10);
    
    size_ expected_size = _mask.count();
    if (parameters.size() != expected_size)
      throw std::range_error("Incorrect item count");
    
    _dat.fill(0);
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  test_case("IXFFDihedral") {
    using FTyp = DihedralFunctionType;
    using TDhd = test::TestDihedralType;
    using AllowEnum = indigox::test::TestDihedralType::AllowEnum;
    using StoreEnum = indigox::test::TestDihedralType::StoreEnum;
    
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
// == IXFFAngleType Serialisation =============================================
// ============================================================================
  
  template <class Archive>
  void IXFFAngleType::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("mask", _mask));
  }
  
  template <class Archive>
  void IXFFAngleType::load_and_construct(Archive &archive,
                                            cereal::construct<IXFFAngleType> &construct,
                                            const uint32_t) {
    construct(AngleFunctionType::Empty, 0, std::initializer_list<float_>());
    archive(INDIGOX_SERIAL_NVP("type", construct->_type),
            INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("mask", construct->_mask));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFAngleType);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAngle serilisation", T, ixffang_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestAngleType saved(AngleFunctionType::Harmonic, 2, {0.1,0.2});
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestAngleType loaded(AngleFunctionType::Empty, 0, {});
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
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffang_serial, ixserial<IXFFAngleType>);
  
// ============================================================================
// == IXFFAngleType Construction ==============================================
// ============================================================================
  
  IXFFAngleType::IXFFAngleType(AngleFunctionType type, int_ id,
                               std::initializer_list<float_> parameters)
  : _type(type), _id(id) {
    // Allow force constant and ideal angle
    if (_type == AngleFunctionType::Harmonic) _mask.from_uint32(3);
    // Allow force constant and ideal angle
    if (_type == AngleFunctionType::CosineHarmonic) _mask.from_uint32(3);
    
    size_ expected_size = _mask.count();
    if (parameters.size() != expected_size)
      throw std::range_error("Incorrect item count");
    
    _dat.fill(0);
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  test_case("IXFFAngle") {
    using FTyp = AngleFunctionType;
    using TAng = test::TestAngleType;
    using AllowEnum = indigox::test::TestAngleType::AllowEnum;
    using StoreEnum = indigox::test::TestAngleType::StoreEnum;
    
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
// == IXFFBondType Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXFFBondType::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("data", _dat),
            INDIGOX_SERIAL_NVP("mask", _mask));
  }
  
  template <class Archive>
  void IXFFBondType::load_and_construct(Archive &archive,
                                        cereal::construct<IXFFBondType> &construct,
                                        const uint32_t) {
    construct(BondFunctionType::Empty, 0, std::initializer_list<float_>());
    archive(INDIGOX_SERIAL_NVP("type", construct->_type),
            INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("data", construct->_dat),
            INDIGOX_SERIAL_NVP("mask", construct->_mask));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFBondType);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFBond serilisation", T, ixffbnd_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestBondType saved(BondFunctionType::Harmonic, 2, {0.1,0.2});
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestBondType loaded(BondFunctionType::Empty, 0, {});
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
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffbnd_serial, ixserial<IXFFBondType>);
  
// ============================================================================
// == IXFFBondType Construction ===============================================
// ============================================================================
  
  IXFFBondType::IXFFBondType(BondFunctionType type, int_ id,
                             std::initializer_list<float_> parameters)
  : _type(type), _id(id) {
    // Allow force constant and ideal length
    if (_type == BondFunctionType::Harmonic) _mask.from_uint32(3);
    // Allow force constant and ideal length
    if (_type == BondFunctionType::Quartic) _mask.from_uint32(3);
    
    size_ expected_size = _mask.count();
    if (parameters.size() != expected_size)
      throw std::range_error("Incorrect item count");
    
    _dat.fill(0);
    std::copy(parameters.begin(), parameters.end(), _dat.begin());
  }
  
  test_case("IXFFBond") {
    using FTyp = BondFunctionType;
    using TBnd = test::TestBondType;
    using AllowEnum = indigox::test::TestBondType::AllowEnum;
    using StoreEnum = indigox::test::TestBondType::StoreEnum;
    
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
// == IXFFAtomType Serialisation ==============================================
// ============================================================================
  
  template <class Archive>
  void IXFFAtomType::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("id", _id),
            INDIGOX_SERIAL_NVP("name", _name));
  }
  
  template <class Archive>
  void IXFFAtomType::load_and_construct(Archive &archive,
                                        cereal::construct<IXFFAtomType> &construct,
                                        const uint32_t) {
    construct(0, "blank");
    archive(INDIGOX_SERIAL_NVP("id", construct->_id),
            INDIGOX_SERIAL_NVP("name", construct->_name));
  }
  INDIGOX_SERIALISE_SPLIT(IXFFAtomType);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAtom serialisation", T, ixffatm_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestAtomType saved(12, "tester");
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp));
    }
    
    test::TestAtomType loaded(0,"");
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded.imp));
    }
    
    check_eq(saved.get_id(), loaded.get_id());
    check_eq(saved.get_name(), loaded.get_name());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffatm_serial, ixserial<IXFFAtomType>);
  
// ============================================================================
// == IXFFAtomType Construction ===============================================
// ============================================================================
  
  IXFFAtomType::IXFFAtomType(int_ id, string_ name) : _id(id), _name(name) { }
  
  test_case("IXFFAtom") {
    subcase("Construction check") {
      check_nothrow(test::TestAtomType t(0, ""));
      
      test::TestAtomType testatm(2, "CH0");
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
    construct(ForcefieldFamily::Empty, "");
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
  
  IXForcefield::IXForcefield(ForcefieldFamily family, string_ name)
  : _family(family), _name(name) {
    
    if (family == ForcefieldFamily::GROMOS) {
      _bnds.emplace(BondFunctionType::Harmonic, BondTypes::mapped_type());
      _bnds.emplace(BondFunctionType::Quartic, BondTypes::mapped_type());
      _angs.emplace(AngleFunctionType::Harmonic, AngleTypes::mapped_type());
      _angs.emplace(AngleFunctionType::CosineHarmonic, AngleTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::Proper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::Improper, DihedralTypes::mapped_type());
    }
    if (family == ForcefieldFamily::Other) {
      _bnds.emplace(BondFunctionType::Harmonic, BondTypes::mapped_type());
      _bnds.emplace(BondFunctionType::Quartic, BondTypes::mapped_type());
      _bnds.emplace(BondFunctionType::Morse, BondTypes::mapped_type());
      _bnds.emplace(BondFunctionType::Cubic, BondTypes::mapped_type());
      _bnds.emplace(BondFunctionType::FENE, BondTypes::mapped_type());
      _angs.emplace(AngleFunctionType::Harmonic, AngleTypes::mapped_type());
      _angs.emplace(AngleFunctionType::CosineHarmonic, AngleTypes::mapped_type());
      _angs.emplace(AngleFunctionType::UreyBradley, AngleTypes::mapped_type());
      _angs.emplace(AngleFunctionType::Quartic, AngleTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::Proper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::Improper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::RyckaertBellemans, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::PeriodicImproper, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::Fourier, DihedralTypes::mapped_type());
      _dhds.emplace(DihedralFunctionType::Restricted, DihedralTypes::mapped_type());
    }
  }
  
// ============================================================================
// == IXForcefield Adding types ===============================================
// ============================================================================
  
  FFAtomType IXForcefield::NewAtomType(int_ id, string_ name) {
    if (GetAtomType(id) || GetAtomType(name)) return FFAtomType();
    _atms.emplace_back(new IXFFAtomType(id, name));
    return _atms.back();
  }
  
  FFBondType IXForcefield::NewBondType(BondFunctionType type, int_ id,
                                       std::initializer_list<float_> parameters) {
    if (_bnds.find(type) == _bnds.end()) return FFBondType();
    if (GetBondType(type, id)) return FFBondType();
    _bnds[type].emplace_back(new IXFFBondType(type, id, parameters));
    return _bnds[type].back();
  }
  
  FFAngleType IXForcefield::NewAngleType(AngleFunctionType type, int_ id,
                                       std::initializer_list<float_> parameters) {
    if (_angs.find(type) == _angs.end()) return FFAngleType();
    if (GetAngleType(type, id)) return FFAngleType();
    _angs[type].emplace_back(new IXFFAngleType(type, id, parameters));
    return _angs[type].back();
  }
  
  FFDihedralType IXForcefield::NewDihedralType(DihedralFunctionType type, int_ id,
                                       std::initializer_list<float_> parameters) {
    if (_dhds.find(type) == _dhds.end()) return FFDihedralType();
    if (GetDihedralType(type, id)) return FFDihedralType();
    _dhds[type].emplace_back(new IXFFDihedralType(type, id, parameters));
    return _dhds[type].back();
  }
  
// ============================================================================
// == IXForcefield Getting types ==============================================
// ============================================================================
  
  FFAtomType IXForcefield::GetAtomType(string_ name) const {
    auto pred = [&name](const FFAtomType& t) { return t->GetName() == name; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? FFAtomType() : *pos;
  }
  
  FFAtomType IXForcefield::GetAtomType(int_ id) const {
    auto pred = [&id](const FFAtomType& t) { return t->GetID() == id; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? FFAtomType() : *pos;
  }
  
  FFBondType IXForcefield::GetBondType(BondFunctionType type, int_ id) const {
    if (_bnds.find(type) == _bnds.end()) return FFBondType();
    auto pred = [&id](const FFBondType& t) { return t->GetID() == id; };
    auto pos = std::find_if(_bnds.at(type).begin(), _bnds.at(type).end(), pred);
    return (pos == _bnds.at(type).end()) ? FFBondType() : *pos;
  }
  
  FFAngleType IXForcefield::GetAngleType(AngleFunctionType type, int_ id) const {
    if (_angs.find(type) == _angs.end()) return FFAngleType();
    auto pred = [&id](const FFAngleType& t) { return t->GetID() == id; };
    auto pos = std::find_if(_angs.at(type).begin(), _angs.at(type).end(), pred);
    return (pos == _angs.at(type).end()) ? FFAngleType() : *pos;
  }
  
  FFDihedralType IXForcefield::GetDihedralType(DihedralFunctionType type, int_ id) const {
    if (_dhds.find(type) == _dhds.end()) return FFDihedralType();
    auto pred = [&id](const FFDihedralType& t) { return t->GetID() == id; };
    auto pos = std::find_if(_dhds.at(type).begin(), _dhds.at(type).end(), pred);
    return (pos == _dhds.at(type).end()) ? FFDihedralType() : *pos;
  }
  
  test_case("IXForcefield") {
    test::TestForcefield ff(ForcefieldFamily::GROMOS, "Test");
    ff.imp = GenerateGROMOS54A7();
    check_eq(0, ff.NumAtomTypes());
    check_eq(0, ff.NumBondTypes());
    check_eq(0, ff.NumAngleTypes());
    check_eq(0, ff.NumDihedralTypes());
  }
  
  Forcefield GenerateGROMOS54A7() {
    Forcefield ff = std::make_shared<IXForcefield>(ForcefieldFamily::GROMOS,
                                                   "GROMOS 54A7");
    // Add bond types
    ff->NewQuarticBondType(1, 15.7e6, 0.1);
    ff->NewHarmonicBondType(1, 0.314e6, 0.1);
    ff->NewQuarticBondType(40, 8.12e6, 0.1758);
    ff->NewHarmonicBondType(40, 0.502e6, 0.1758);
    
    // Add Angle types
    ff->NewCosineHarmonicAngleType(2, 420, 90.0);
    ff->NewHarmonicAngleType(2, 0.128, 90.0);
    ff->NewCosineHarmonicAngleType(31, 700, 122.0);
    ff->NewCosineHarmonicAngleType(31, 0.153, 122.0);
    
    // Add Improper types
    ff->NewImproperDihedralType(2, 0.102, 35.36439);
    ff->NewImproperDihedralType(4, 0.0510, 180.0);
    
    // Add Proper types
    ff->NewProperDihedralType(30, 3.9, 0, 3);
    ff->NewProperDihedralType(15, 41.8, 180, 2);
    
    // Add Ato types
    ff->NewAtomType(37, "NA+");
    ff->NewAtomType(44, "ODmso");
    return ff;
  }
  
  test_suite_close();
}
