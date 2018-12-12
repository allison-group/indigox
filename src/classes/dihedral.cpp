#include <indigox/classes/atom.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/utils/serialise.hpp>

#include <memory>
#include <vector>

namespace indigox {

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid dihedral instance")
#else
#define _sanity_check_(x)
#endif

  test_suite_open("Dihedral");

  // =======================================================================
  // == SERIALISATION ======================================================
  // =======================================================================

  template <typename Archive>
  void Dihedral::Impl::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("molecule", molecule),
            INDIGOX_SERIAL_NVP("tag", tag),
            INDIGOX_SERIAL_NVP("unique_id", unique_id),
            INDIGOX_SERIAL_NVP("types", forcefield_types),
            INDIGOX_SERIAL_NVP("priority", priority));
  }

  template <typename Archive>
  void Dihedral::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }

  // =======================================================================
  // == OPERATORS ==========================================================
  // =======================================================================

  bool Dihedral::operator==(const Dihedral &dhd) const {
    return m_data == dhd.m_data;
  }

  bool Dihedral::operator<(const Dihedral &dhd) const {
    _sanity_check_(*this);
    _sanity_check_(dhd);
    return m_data->unique_id < dhd.m_data->unique_id;
  }

  bool Dihedral::operator>(const Dihedral &dhd) const {
    _sanity_check_(*this);
    _sanity_check_(dhd);
    return m_data->unique_id > dhd.m_data->unique_id;
  }

  // =======================================================================
  // == CONSTRUCTION =======================================================
  // =======================================================================

  Dihedral::Impl::Impl(const Atom &a, const Atom &b, const Atom &c,
                       const Atom &d, const Molecule &mol)
      : atoms({a, b, c, d}), molecule(mol), priority(0) {
  }

  Dihedral::Dihedral(const Atom &a, const Atom &b, const Atom &c, const Atom &d,
                     const Molecule &mol)
      : m_data(std::make_shared<Impl>(a, b, c, d, mol)) {
  }

  void Dihedral::Reset() {
    m_data->atoms.fill(Atom());
    m_data->molecule = Molecule();
    m_data->tag = INT64_MIN;
    m_data->unique_id = INT64_MIN;
    m_data->forcefield_types.clear();
    m_data->priority = INT32_MIN;
  }

  // =======================================================================
  // == STATE CHECKING =====================================================
  // =======================================================================

  int64_t Dihedral::NumTypes() const {
    _sanity_check_(*this);
    return (int64_t)(m_data->forcefield_types.size());
  }

  bool Dihedral::HasType() const {
    _sanity_check_(*this);
    return !m_data->forcefield_types.empty();
  }

  // =======================================================================
  // == STATE GETTING ======================================================
  // =======================================================================

  int64_t Dihedral::GetTag() const {
    _sanity_check_(*this);
    return m_data->tag;
  }

  int64_t Dihedral::GetID() const {
    _sanity_check_(*this);
    return m_data->unique_id;
  }

  const Molecule &Dihedral::GetMolecule() const {
    _sanity_check_(*this);
    return m_data->molecule;
  }

  const Dihedral::DihedralAtoms &Dihedral::GetAtoms() const {
    _sanity_check_(*this);
    return m_data->atoms;
  }

  int64_t Dihedral::GetIndex() const {
    _sanity_check_(*this);
    if (!m_data->molecule) {
      return -1;
    }
    auto mol_dihedrals = m_data->molecule.GetDihedrals();
    auto idx = std::find(mol_dihedrals.begin(), mol_dihedrals.end(), *this);
    return (idx == mol_dihedrals.end())
               ? -1
               : std::distance(mol_dihedrals.begin(), idx);
  }

  const Dihedral::DihedralTypes &Dihedral::GetTypes() const {
    _sanity_check_(*this);
    return m_data->forcefield_types;
  }

  int32_t Dihedral::GetPriority() const {
    _sanity_check_(*this);
    return m_data->priority;
  }

  // =======================================================================
  // == STATE SETTING ======================================================
  // =======================================================================

  void Dihedral::SetTag(int64_t tag) {
    _sanity_check_(*this);
    m_data->tag = tag;
  }

  void Dihedral::SetTypes(const DihedralTypes &types) {
    _sanity_check_(*this);
    m_data->forcefield_types.assign(types.begin(), types.end());
  }

  void Dihedral::AddType(const FFDihedral &type) {
    _sanity_check_(*this);
    if (m_data->molecule) {
      if (!m_data->molecule.HasForcefield()) {
        m_data->molecule.SetForcefield(type.GetForcefield());
      } else if (m_data->molecule.GetForcefield() != type.GetForcefield()) {
        throw std::runtime_error(
            "Provided dihedral type does not match molecule's forcefield");
      }
    }

    m_data->forcefield_types.emplace_back(type);
    std::sort(m_data->forcefield_types.begin(), m_data->forcefield_types.end());
  }

  void Dihedral::RemoveType(const FFDihedral &type) {
    _sanity_check_(*this);
    auto pos = std::find(m_data->forcefield_types.begin(),
                         m_data->forcefield_types.end(), type);
    if (pos != m_data->forcefield_types.end())
      m_data->forcefield_types.erase(pos);
  }

  // =======================================================================
  // == OUTPUTTING =========================================================
  // =======================================================================

  /*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXDihedral serialisation", T,
    ixdihed_serial) { using In = typename T::t1; using Out = typename
    cereal::traits::detail::get_output_from_input<In>::type;
      test::DihedralTestFixture fixture;
      Dihedral saved = fixture.dhd.imp;

      saved->SetTag(45);
      saved->SetType(fixture.fftype);
      std::ostringstream os;
      {
        Out oar(os);
        check_nothrow(oar(saved, saved->GetAtoms(), saved->GetMolecule(),
                          fixture.fftype));
      }

      Dihedral loaded;
      test::TestDihedral::Atoms atoms;
      Molecule mol;
      FFDihedral ff_loaded;
      std::istringstream is(os.str());
      {
        In iar(is);
        check_nothrow(iar(loaded, atoms, mol, ff_loaded));
      }
      fixture.dhd.imp = loaded;
      check(!fixture.dhd.get_atms()[0].expired());
      check(!fixture.dhd.get_atms()[1].expired());
      check(!fixture.dhd.get_atms()[2].expired());
      check(!fixture.dhd.get_atms()[3].expired());
      check(!fixture.dhd.get_mol().expired());
      check_eq(saved->GetTag(), loaded->GetTag());
      check_eq(atoms.first->GetTag(), saved->GetAtoms().first->GetTag());
      check_eq(atoms.second->GetTag(), saved->GetAtoms().second->GetTag());
      check_eq(atoms.third->GetTag(), saved->GetAtoms().third->GetTag());
      check_eq(atoms.fourth->GetTag(), saved->GetAtoms().fourth->GetTag());
      check_eq(atoms, loaded->GetAtoms());
      check_eq(ff_loaded, loaded->GetType());
    }
    DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixdihed_serial,
    ixserial<IXDihedral>);
   */

  /*  test_case_fixture(test::DihedralTestFixture, "IXDihedral construction") {
      check_nothrow(test::TestDihedral tdhd(a,b,c,d,mol));
      test::TestDihedral tdhd(a,b,c,d,mol);
      check_eq(mol, tdhd.get_mol().lock());
      check_eq(0, tdhd.get_tag());
      check_eq(a, tdhd.get_atms()[0].lock());
      check_eq(b, tdhd.get_atms()[1].lock());
      check_eq(c, tdhd.get_atms()[2].lock());
      check_eq(d, tdhd.get_atms()[3].lock());
      check_eq(FFDihedral(), tdhd.get_type());

      // Check unique IDs correctly update
      test::TestDihedral dhd1(a,b,c,d,mol);
      test::TestDihedral dhd2(a,b,c,d,mol);
      check_ne(dhd1.GetUniqueID(), dhd2.GetUniqueID());
      check_eq(dhd1.GetUniqueID() + 1, dhd2.GetUniqueID());
    }
   */

  /*  test_case_fixture(test::DihedralTestFixture, "IXDihedral getting and
    setting") {
      // check no throwing
      check_nothrow(dhd.GetTag());
      check_nothrow(dhd.SetTag(72));
      check_nothrow(dhd.GetMolecule());
      check_nothrow(dhd.GetAtoms());
      check_nothrow(dhd.SwapOrder());
      check_nothrow(dhd.NumAtoms());
      check_nothrow(dhd.SetType(fftype));

      // Check correctness of gets
      check_eq(72, dhd.GetTag());
      check_eq(72, dhd.get_tag());
      check_eq(mol, dhd.GetMolecule());
      check_eq(mol, dhd.get_mol().lock());
      check_eq(stdx::make_quad(d, c, b, a), dhd.GetAtoms());
      check_eq(dhd.GetTag(), dhd.GetIndex());
      check_eq(fftype, dhd.GetType());

      // Check correct no owning of molecule
      mol.reset();
      check(dhd.get_mol().expired());
      check_eq(Molecule(), dhd.GetMolecule());

      // Check still gets tag when mol is dead
      check_eq(dhd.GetTag(), dhd.GetIndex());

      // Check num atoms doesn't depend on being active
      check_eq(4, dhd.NumAtoms());
      a.reset(); b.reset(); c.reset(); d.reset();
      check_eq(4, dhd.NumAtoms());
    }
   */

  std::ostream &operator<<(std::ostream &os, const Dihedral &dhd) {
    _sanity_check_(dhd);
    os << "Dihedral(" << dhd.m_data->atoms[0] + 1 << ", "
       << dhd.m_data->atoms[1] + 1 << ", " << dhd.m_data->atoms[2] + 1 << ", "
       << dhd.m_data->atoms[3] + 1 << ")";
    return os;
  }

  /*  test_case_fixture(test::DihedralTestFixture, "IXDihedral printing
    methods") {
      // Check ordering is correct
      std::stringstream ss; ss << dhd.imp;
      check_eq("Dihedral(Atom(0, C), Atom(1, O), Atom(2, F), Atom(3, N))",
               dhd.ToString());
      check_eq("Dihedral(0, 1, 2, 3)", ss.str());
      dhd.SwapOrder();
      ss.str(""); ss << dhd.imp;
      check_eq("Dihedral(Atom(3, N), Atom(2, F), Atom(1, O), Atom(0, C))",
               dhd.ToString());
      check_eq("Dihedral(3, 2, 1, 0)", ss.str());

      // Check having bad atoms is handled correctly
      a.reset(); d.reset();
      ss.str(""); ss << dhd.imp;
      check_eq("Dihedral(MALFORMED)", dhd.ToString());
      check_eq("Dihedral(, 2, 1, )", ss.str());

      // Check empty dihedral to ostream does nothing
      ss.str("");
      check_nothrow(ss << Dihedral());
      check_eq("", ss.str());
    }
   */

  test_suite_close();
} // namespace indigox
