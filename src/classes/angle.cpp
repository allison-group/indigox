#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/doctest_proxy.hpp>

#include <array>
#include <memory>

namespace indigox {
  test_suite_open("Angle");

  // =======================================================================
  // == IMPLEMENTATION =====================================================
  // =======================================================================

  struct Angle::Impl {
    AngleAtoms atoms;
    Molecule molecule;
    int64_t tag;
    int64_t unique_id;
    FFAngle forcefield_type;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t);
  };

  // =======================================================================
  // == SERIALISATION ======================================================
  // =======================================================================

  template <typename Archive>
  void Angle::Impl::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("molecule", molecule),
            INDIGOX_SERIAL_NVP("tag", tag),
            INDIGOX_SERIAL_NVP("unique_id", unique_id),
            INDIGOX_SERIAL_NVP("type", forcefield_type));
  }

  template <typename Archive>
  void Angle::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(Angle);

  // =======================================================================
  // == OPERATORS ==========================================================
  // =======================================================================

  bool Angle::operator==(const Angle &ang) const {
    if (ang.m_data == m_data) {
      return true;
    } else {
      return false;
    }

    // Do a more indepth comparison if they're not the exact same angle
    if (ang.m_data->forcefield_type != m_data->forcefield_type) {
      return false;
    }
    if (ang.m_data->atoms[1] != m_data->atoms[1]) {
      return false;
    }

    bool fwd_comp = (ang.m_data->atoms[0] == m_data->atoms[0]) &&
                    (ang.m_data->atoms[2] == m_data->atoms[2]);
    bool bwd_comp = (ang.m_data->atoms[2] == m_data->atoms[0]) &&
                    (ang.m_data->atoms[0] == m_data->atoms[2]);
    return (fwd_comp || bwd_comp);
  }

  bool Angle::operator<(const Angle &ang) const {
    return m_data->unique_id < ang.m_data->unique_id;
  }
  bool Angle::operator>(const Angle &ang) const {
    return m_data->unique_id > ang.m_data->unique_id;
  }

  // =======================================================================
  // == CONSTRUCTION =======================================================
  // =======================================================================

  Angle::Angle(const Atom &a, const Atom &b, const Atom &c, const Molecule &m)
      : m_data(std::make_shared<Impl>()) {
    m_data->atoms[0] = a;
    m_data->atoms[1] = b;
    m_data->atoms[2] = c;
    m_data->molecule = m;
  }

  // =======================================================================
  // == STATE CHECKING =====================================================
  // =======================================================================

  bool Angle::HasType() const {
    return bool(*this) ? bool(m_data->forcefield_type) : false;
  }

  // =======================================================================
  // == STATE GETTING ======================================================
  // =======================================================================

  int64_t Angle::GetTag() const {
    return *this ? m_data->tag
                 : throw std::runtime_error(
                       "Attempting to access data from invalid angle instance");
  }

  int64_t Angle::GetID() const {
    return *this ? m_data->unique_id
                 : throw std::runtime_error(
                       "Attempting to access data from invalid angle instance");
  }

  const Molecule &Angle::GetMolecule() const {
    return *this ? m_data->molecule
                 : throw std::runtime_error(
                       "Attempting to access data from invalid angle instance");
  }

  const Angle::AngleAtoms &Angle::GetAtoms() const {
    return *this ? m_data->atoms
                 : throw std::runtime_error(
                       "Attempting to access data from invalid angle instance");
  }

  const FFAngle &Angle::GetType() const {
    return *this ? m_data->forcefield_type
                 : throw std::runtime_error(
                       "Attempting to access data from invalid angle instance");
  }

  int64_t Angle::GetIndex() const {
    if (!(*this)) {
      throw std::runtime_error(
          "Attempting to access data from invalid angle instance");
    }
    if (!m_data->molecule) {
      return -1;
    }
    auto mol_angles = m_data->molecule.GetAngles();
    auto idx = std::find(mol_angles.begin(), mol_angles.end(), *this);
    return (idx == mol_angles.end()) ? -1
                                     : std::distance(mol_angles.begin(), idx);
  }

  // =======================================================================
  // == STATE SETTING ======================================================
  // =======================================================================

  void Angle::SetType(const FFAngle &type) {
    if (!(*this)) {
      throw std::runtime_error(
          "Attempting to access data from invalid angle instance");
    }
    if (m_data->molecule) {
      if (!m_data->molecule.HasForcefield()) {
        m_data->molecule.SetForcefield(type.GetForcefield());
      } else if (m_data->molecule.GetForcefield() != type.GetForcefield()) {
        throw std::runtime_error(
            "Provided angle type does not match molecule's forcefield");
      }
    }
    m_data->forcefield_type = type;
  }

  // =======================================================================
  // == OUTPUTTING =========================================================
  // =======================================================================

  std::ostream &operator<<(std::ostream &os, const Angle &ang) {
    if (ang) {
      const Angle::AngleAtoms &atoms = ang.GetAtoms();
      os << "Angle(" << atoms[0].GetIndex() + 1 << ", "
         << atoms[1].GetIndex() + 1 << ", " << atoms[2].GetIndex() + 1 << ")";
    }
    return os;
  }

  /*    DOCTEST_TEST_CASE_TEMPLATE_DEFINE("Angle serialisation", T,
   ixangle_serial) { using In = typename T::t1; using Out = typename
   cereal::traits::detail::get_output_from_input<In>::type; FFAngle ffang =
   test::CreateGenericTestFFAngle().imp; test::AngleTestFixture fixture; Angle
   saved = fixture.ang.imp;

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

  /*    test_case_fixture(test::AngleTestFixture, "IXAngle getting and setting")
     { FFAngle ffang = test::CreateGenericTestFFAngle().imp;

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

  test_suite_close();
} // namespace indigox
