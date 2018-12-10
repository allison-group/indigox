#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/utils/serialise.hpp>

#include <array>
#include <memory>

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid bond instance")
#else
#define _sanity_check_(x)
#endif

namespace indigox {

  test_suite_open("Bond");

  // =======================================================================
  // == IMPLEMENTATION =====================================================
  // =======================================================================

  struct Bond::Impl {
    BondAtoms atoms;
    Molecule molecule;
    int64_t tag;
    int64_t unique_id;
    BondOrder order;
    BondStereo stereochemistry;
    FFBond forcefield_type;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

    Impl() = default;
    Impl(const Atom &a, const Atom &b, const Molecule &mol, BondOrder o);
  };

  // =======================================================================
  // == SERIALISATION ======================================================
  // =======================================================================

  template <typename Archive>
  void Bond::Impl::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("molecule", molecule),
            INDIGOX_SERIAL_NVP("tag", tag),
            INDIGOX_SERIAL_NVP("unique_id", unique_id),
            INDIGOX_SERIAL_NVP("order", order),
            INDIGOX_SERIAL_NVP("sterochemistry", stereochemistry),
            INDIGOX_SERIAL_NVP("type", forcefield_type));
  }

  template <typename Archive>
  void Bond::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(Bond);

  // =======================================================================
  // == OPERATORS ==========================================================
  // =======================================================================

  bool Bond::operator==(const Bond &bnd) const {
    return m_data == bnd.m_data;
  }

  bool Bond::operator<(const Bond &bnd) const {
    _sanity_check_(*this);
    _sanity_check_(bnd);
    return m_data->unique_id < bnd.m_data->unique_id;
  }

  bool Bond::operator>(const Bond &bnd) const {
    _sanity_check_(*this);
    _sanity_check_(bnd);
    return m_data->unique_id > bnd.m_data->unique_id;
  }

  // =======================================================================
  // == CONSTRUCTION =======================================================
  // =======================================================================

  Bond::Impl::Impl(const Atom &a, const Atom &b, const Molecule &mol,
                   BondOrder o)
      : atoms({a, b}), molecule(mol), order(o),
        stereochemistry(BondStereo::UNDEFINED) {
  }

  Bond::Bond(const Atom &a, const Atom &b, const Molecule &m, BondOrder o)
      : m_data(std::make_shared<Impl>(a, b, m, o)) {
  }

  // =======================================================================
  // == STATE CHECKING =====================================================
  // =======================================================================

  bool Bond::HasType() const {
    _sanity_check_(*this);
    return bool(m_data->forcefield_type);
  }

  // =======================================================================
  // == STATE GETTING ======================================================
  // =======================================================================

  int64_t Bond::GetTag() const {
    _sanity_check_(*this);
    return m_data->tag;
  }

  int64_t Bond::GetID() const {
    _sanity_check_(*this);
    return m_data->unique_id;
  }

  const Molecule &Bond::GetMolecule() const {
    _sanity_check_(*this);
    return m_data->molecule;
  }

  BondOrder Bond::GetOrder() const {
    _sanity_check_(*this);
    return m_data->order;
  }

  BondStereo Bond::GetStereochemistry() const {
    _sanity_check_(*this);
    return m_data->stereochemistry;
  }

  const Bond::BondAtoms &Bond::GetAtoms() const {
    _sanity_check_(*this);
    return m_data->atoms;
  }

  int64_t Bond::GetIndex() const {
    _sanity_check_(*this);
    if (!m_data->molecule) {
      return -1;
    }
    auto mol_bonds = m_data->molecule.GetBonds();
    auto idx = std::find(mol_bonds.begin(), mol_bonds.end(), *this);
    return (idx == mol_bonds.end()) ? -1
                                    : std::distance(mol_bonds.begin(), idx);
  }

  const FFBond &Bond::GetType() const {
    _sanity_check_(*this);
    return m_data->forcefield_type;
  }

  // =======================================================================
  // == STATE SETTING ======================================================
  // =======================================================================

  void Bond::SetTag(int64_t tag) {
    _sanity_check_(*this);
    m_data->tag = tag;
  }

  void Bond::SetOrder(BondOrder order) {
    _sanity_check_(*this);
    m_data->order = order;
  }

  void Bond::SetStereochemistry(BondStereo stereo) {
    _sanity_check_(*this);
    m_data->stereochemistry = stereo;
  }

  void Bond::SetType(const FFBond &type) {
    _sanity_check_(*this);
    if (m_data->molecule) {
      if (!m_data->molecule.HasForcefield()) {
        m_data->molecule.SetForcefield(type.GetForcefield());
      } else if (m_data->molecule.GetForcefield() != type.GetForcefield()) {
        throw std::runtime_error(
            "Provided bond type does not match molecule's forcefield");
      }
    }
    m_data->forcefield_type = type;
  }

  // =======================================================================
  // == OUTPUTTING =========================================================
  // =======================================================================

  std::ostream &operator<<(std::ostream &os, const Bond &bnd) {
    if (bnd) {
      os << "Bond(" << bnd.m_data->atoms[0].GetIndex() + 1 << ", "
         << bnd.m_data->atoms[1].GetIndex() + 1 << ")";
    }
    return os;
  }

  /*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXBond serialisation", T,
    ixbond_serial) { using In = typename T::t1; using Out = typename
    cereal::traits::detail::get_output_from_input<In>::type;
      test::BondTestFixture fixture;
      Bond saved = fixture.bnd.imp;

      saved->SetTag(12);
      saved->SetOrder(BondOrder::TRIPLE);
      saved->SetAromaticity(true);
      saved->SetStereochemistry(BondStereo::Z);
      saved->SetType(fixture.fftype);

      std::ostringstream os;
      {
        Out oar(os);
        // atoms and molecule saved so they can be kept alive on load
        check_nothrow(oar(saved, saved->GetSourceAtom(), saved->GetTargetAtom(),
                          saved->GetMolecule(), fixture.fftype));
      }

      Bond loaded;
      Atom source, target;
      Molecule mol;
      FFBond typ_loaded;
      std::istringstream is(os.str());
      {
        In iar(is);
        check_nothrow(iar(loaded, source, target, mol, typ_loaded));
      }
      fixture.bnd.imp = loaded;
      check(!fixture.bnd.get_atms()[0].expired());
      check(!fixture.bnd.get_atms()[1].expired());
      check(!fixture.bnd.get_mol().expired());
      check_eq(saved->GetTag(), loaded->GetTag());
      check_eq(saved->GetOrder(), loaded->GetOrder());
      check_eq(saved->GetAromaticity(), loaded->GetAromaticity());
      check_eq(saved->GetStereochemistry(), loaded->GetStereochemistry());
      check_eq(source->GetTag(), saved->GetSourceAtom()->GetTag());
      check_eq(target->GetTag(), saved->GetTargetAtom()->GetTag());
      check_eq(source, loaded->GetSourceAtom());
      check_eq(target, loaded->GetTargetAtom());
      check_eq(mol, loaded->GetMolecule());
      check_eq(typ_loaded, loaded->GetType());
    }
    DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixbond_serial, ixserial<IXBond>);
   */

  /*  test_case_fixture(test::BondTestFixture, "IXBond construction") {
      check_nothrow(test::TestBond tbnd(a,b,mol));
      test::TestBond tbnd(a,b,mol);
      check_eq(mol, tbnd.get_mol().lock());
      check_eq(0, tbnd.get_tag());
      check_eq(BondOrder::UNDEFINED, tbnd.get_order());
      check_eq(false, tbnd.get_aromatic());
      check_eq(BondStereo::UNDEFINED, tbnd.get_stereo());
      check_eq(a, tbnd.get_atms()[0].lock());
      check_eq(b, tbnd.get_atms()[1].lock());
      check_eq(FFBond(), tbnd.get_type());

      // Check unique IDs correctly update
      test::TestBond bnd1(a,b,mol);
      test::TestBond bnd2(a,b,mol);
      check_ne(bnd1.GetUniqueID(), bnd2.GetUniqueID());
      check_eq(bnd1.GetUniqueID() + 1, bnd2.GetUniqueID());
    }
   */

  /*  test_case_fixture(test::BondTestFixture, "IXBond getting and setting") {
      // Check no throwing
      check_nothrow(bnd.SetTag(12));
      check_nothrow(bnd.SetOrder(BondOrder::SINGLE));
      check_nothrow(bnd.SetAromaticity(true));
      check_nothrow(bnd.SetStereochemistry(BondStereo::NONE));
      check_nothrow(bnd.SwapSourceTarget());
      check_nothrow(bnd.SetType(fftype));

      // Check correctness of gets
      check_eq(12, bnd.GetTag());
      check_eq(12, bnd.get_tag());
      check_eq(mol, bnd.GetMolecule());
      check_eq(mol, bnd.get_mol().lock());
      check_eq(BondOrder::SINGLE, bnd.GetOrder());
      check_eq(BondOrder::SINGLE, bnd.get_order());
      check_eq(b, bnd.GetSourceAtom());
      check_eq(a, bnd.GetTargetAtom());
      check_eq(a, bnd.get_atms()[1].lock());
      check_eq(b, bnd.get_atms()[0].lock());
      check_eq(true, bnd.GetAromaticity());
      check_eq(true, bnd.get_aromatic());
      check_eq(BondStereo::NONE, bnd.GetStereochemistry());
      check_eq(BondStereo::NONE, bnd.get_stereo());
      check_eq(std::make_pair(b,a), bnd.GetAtoms());
      check_eq(bnd.GetTag(), bnd.GetIndex());
      check_eq(fftype, bnd.GetType());

      // Check no owning of molecule
      mol.reset();
      check(bnd.get_mol().expired());
      check_eq(Molecule(), bnd.GetMolecule());
      // Check still gets tag when mol is dead
      check_eq(bnd.GetTag(), bnd.GetIndex());

      // Check num atoms doesn't depend on being active
      check_eq(2, bnd.NumAtoms());
      a.reset(); b.reset();
      check_eq(2, bnd.NumAtoms());
    }
   */

  /*  test_case_fixture(test::BondTestFixture, "IXBond printing methods") {
      // Check ordering is correct
      std::stringstream ss;
      check_nothrow(ss << bnd.imp);
      check_eq("Bond(Atom(0, C), Atom(1, O))", bnd.ToString());
      check_eq("Bond(0, 1)", ss.str());
      bnd.SwapSourceTarget();
      ss.str("");
      check_nothrow(ss << bnd.imp);
      check_eq("Bond(Atom(1, O), Atom(0, C))", bnd.ToString());
      check_eq("Bond(1, 0)", ss.str());

      // Check having bad atoms is handled correctly
      a.reset();
      ss.str("");
      check_nothrow(ss << bnd.imp);
      check_eq("Bond(MALFORMED)", bnd.ToString());
      check_eq("Bond(1, )", ss.str());

      // Check empty bond to ostream does nothing
      ss.str("");
      check_nothrow(ss << Bond());
      check_eq("", ss.str());
    }
   */

  test_suite_close();

} // namespace indigox
