#include <sstream>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/angle_test.hpp>
#include <indigox/test/forcefield_test.hpp>

namespace indigox {
  test_suite_open("Angle");
  
  // Serialisation of Angle
  template <typename Archive>
  void Angle::save(Archive &archive, const uint32_t) const {
    std::vector<sAtom> atoms;
    for (wAtom at : _atms) atoms.emplace_back(at.lock());
    archive(INDIGOX_SERIAL_NVP("molecule", _mol.lock()),
            INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("tag", _tag),
            INDIGOX_SERIAL_NVP("type", _type));
  }
  
  template <typename Archive>
  void Angle::load_and_construct(Archive &archive,
                                 cereal::construct<Angle> &construct,
                                 const uint32_t) {
    sMolecule m;
    std::vector<sAtom> atoms;
    archive(INDIGOX_SERIAL_NVP("molecule", m),
            INDIGOX_SERIAL_NVP("atoms", atoms));
    
    construct(*atoms[0], *atoms[1], *atoms[2], *m);
    
    archive(INDIGOX_SERIAL_NVP("tag", construct->_tag),
            INDIGOX_SERIAL_NVP("type", construct->_type));
  }
  INDIGOX_SERIALISE_CONSTRUCT(Angle);
  
/*    DOCTEST_TEST_CASE_TEMPLATE_DEFINE("Angle serialisation", T, ixangle_serial) {
      using In = typename T::t1;
      using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
      FFAngle ffang = test::CreateGenericTestFFAngle().imp;
      test::AngleTestFixture fixture;
      Angle saved = fixture.ang.imp;
  
      saved->SetTag(23);
      saved->SetType(ffang);
  
      std::ostringstream os;
      {
        Out oar(os);
        check_nothrow(oar(saved, saved->GetAtoms(), saved->GetMolecule(),
                          saved->GetType()));
      }
  
      Angle loaded;
      test::TestAngle::Atoms atoms;
      Molecule mol;
      FFAngle loaded_type;
      std::istringstream is(os.str());
      {
        In iar(is);
        check_nothrow(iar(loaded, atoms, mol, loaded_type));
      }
      fixture.ang.imp = loaded;
      check_false(fixture.ang.get_atms()[0].expired());
      check_false(fixture.ang.get_atms()[1].expired());
      check_false(fixture.ang.get_atms()[2].expired());
      check_false(fixture.ang.get_mol().expired());
      check_false(fixture.ang.get_type().use_count() == 0);
      check_eq(saved->GetTag(), loaded->GetTag());
      check_eq(atoms.first->GetTag(), saved->GetAtoms().first->GetTag());
      check_eq(atoms.second->GetTag(), saved->GetAtoms().second->GetTag());
      check_eq(atoms.third->GetTag(), saved->GetAtoms().third->GetTag());
      check_eq(atoms, loaded->GetAtoms());
    }
    DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixangle_serial, ixserial<IXAngle>);
 */
  
  Angle::Angle(Atom& a, Atom& b, Atom& c, Molecule& m)
  : utils::IXCountableObject<Angle>(), _mol(m.weak_from_this()), _tag(0),
  _atms({a.weak_from_this(), b.weak_from_this(), c.weak_from_this()}) { }
  
/*    test_case_fixture(test::AngleTestFixture, "IXAngle construction") {
      check_nothrow(test::TestAngle tang(a,b,c,mol));
      test::TestAngle tang(a,b,c,mol);
      check_eq(mol, tang.get_mol().lock());
      check_eq(0, tang.get_tag());
      check_eq(a, tang.get_atms()[0].lock());
      check_eq(b, tang.get_atms()[1].lock());
      check_eq(c, tang.get_atms()[2].lock());
      check_eq(FFAngle(), tang.get_type());
  
      // Check unique IDs correctly update
      test::TestAngle ang1(a,b,c,mol);
      test::TestAngle ang2(a,b,c,mol);
      check_ne(ang1.GetUniqueID(), ang2.GetUniqueID());
      check_eq(ang1.GetUniqueID() + 1, ang2.GetUniqueID());
    }
 */
  
  size_t Angle::GetIndex() const {
    if (_mol.expired()) return GetTag();
    Molecule& mol = *_mol.lock();
    auto be = mol.GetAngleIters();
    auto pos = std::find(be.first, be.second, shared_from_this());
    if (pos == be.second) return GetTag();
    return std::distance(be.first, pos);
  }
  
  
/*    test_case_fixture(test::AngleTestFixture, "IXAngle getting and setting") {
      FFAngle ffang = test::CreateGenericTestFFAngle().imp;
  
      // Check no throwing
      check_nothrow(ang.GetTag());
      check_nothrow(ang.SetTag(72));
      check_nothrow(ang.GetMolecule());
      check_nothrow(ang.GetAtoms());
      check_nothrow(ang.SwapOrder());
      check_nothrow(ang.NumAtoms());
      check_nothrow(ang.GetType());
      check_nothrow(ang.SetType(ffang));
  
      // Check correctness of gets
      check_eq(72, ang.GetTag());
      check_eq(72, ang.get_tag());
      check_eq(mol, ang.GetMolecule());
      check_eq(mol, ang.get_mol().lock());
      check_eq(stdx::make_triple(c, b, a), ang.GetAtoms());
      check_eq(ang.GetTag(), ang.GetIndex());
      check_eq(ffang, ang.GetType());
  
      // Check correct no owning of molecule
      mol.reset();
      check(ang.get_mol().expired());
      check_eq(Molecule(), ang.GetMolecule());
  
      // Check get index returns tag when molecule dead
      check_eq(ang.GetTag(), ang.GetIndex());
  
      // Check num atoms doesn't depend on being active
      check_eq(3, ang.NumAtoms());
      a.reset(); b.reset(); c.reset();
      check_eq(3, ang.NumAtoms());
    }
 */
  
  std::string Angle::ToString() const {
    std::stringstream ss;
    
    Atom& a = *_atms[0].lock();
    Atom& b = *_atms[1].lock();
    Atom& c = *_atms[2].lock();
    ss << "Angle(" << a.ToString() << ", " << b.ToString() << ", "
       << c.ToString() << ")";
    return ss.str();
  }
  
  std::ostream& operator<<(std::ostream& os, const Angle& ang) {
    auto s = ang.GetAtoms();
    os << "Angle(" << s.first.GetIndex() << ", " << s.second.GetIndex()
       << ", " << s.third.GetIndex() << ")";
    return os;
  }
  
/*    test_case_fixture(test::AngleTestFixture, "IXAngle printing methods") {
      // Check ordering is correct
      std::stringstream ss; ss << ang.imp;
      check_eq("Angle(Atom(0, C), Atom(1, O), Atom(2, N))", ang.ToString());
      check_eq("Angle(0, 1, 2)", ss.str());
      ang.SwapOrder();
      ss.str(""); ss << ang.imp;
      check_eq("Angle(Atom(2, N), Atom(1, O), Atom(0, C))", ang.ToString());
      check_eq("Angle(2, 1, 0)", ss.str());
  
      // Check having bad atoms is handled correctly
      a.reset();
      ss.str(""); ss << ang.imp;
      check_eq("Angle(MALFORMED)", ang.ToString());
      check_eq("Angle(2, 1, )", ss.str());
  
      // Check empty angle to ostream does nothing
      ss.str("");
      check_nothrow(ss << Angle());
      check_eq("", ss.str());
    }
 */
  
  void Angle::SetType(FFAngle& type) { _type = type.weak_from_this(); }
  
  void Angle::Clear() {
    _mol.reset();
    _type.reset();
    _atms.fill(wAtom());
  }
  
  test_suite_close();
}
