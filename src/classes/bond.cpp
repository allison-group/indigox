#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
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

  void Bond::Reset() {
    m_data->atoms.fill(Atom());
    m_data->molecule = Molecule();
    m_data->tag = INT64_MIN;
    m_data->unique_id = INT64_MIN;
    m_data->order = BondOrder::UNDEFINED;
    m_data->stereochemistry = BondStereo::UNDEFINED;
    m_data->forcefield_type = FFBond();
  }

  // =======================================================================
  // == STATE CHECKING =====================================================
  // =======================================================================

  bool Bond::HasType() const {
    _sanity_check_(*this);
    return bool(m_data->forcefield_type);
  }

  bool Bond::IsAmideBond() const {
    _sanity_check_(*this);
    // Check for being an N-C bond
    Atom n, c;
    if (GetAtoms()[0].GetElement() == "N" &&
        GetAtoms()[1].GetElement() == "C") {
      n = GetAtoms()[0];
      c = GetAtoms()[1];
    }
    if (GetAtoms()[1].GetElement() == "N" &&
        GetAtoms()[0].GetElement() == "C") {
      n = GetAtoms()[1];
      c = GetAtoms()[0];
    }
    if (!n || !c)
      return false;
    if (n.NumBonds() + n.GetImplicitCount() != 3)
      return false;

    // Check for C being O bonded.
    for (Bond bnd : c.GetBonds()) {
      if (bnd.IsCarbonylBond())
        return true;
    }
    return false;
  }

  bool Bond::IsCarbonylBond() const {
    _sanity_check_(*this);
    // Check for being C-O bond
    Atom c, o;
    if (GetAtoms()[0].GetElement() == "O" &&
        GetAtoms()[1].GetElement() == "C") {
      o = GetAtoms()[0];
      c = GetAtoms()[1];
    }
    if (GetAtoms()[1].GetElement() == "O" &&
        GetAtoms()[0].GetElement() == "C") {
      o = GetAtoms()[1];
      c = GetAtoms()[0];
    }
    if (!o || !c)
      return false;
    if (o.NumBonds() + o.GetImplicitCount() != 1)
      return false;
    if (c.NumBonds() + c.GetImplicitCount() != 3)
      return false;
    return true;
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

} // namespace indigox
