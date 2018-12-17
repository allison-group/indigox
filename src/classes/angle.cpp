#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
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
      "Attempting to access data from invalid angle instance")
#else
#define _sanity_check_(x)
#endif

namespace indigox {

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
    return (ang.m_data == m_data);

    /* // Do a more indepth comparison if they're not the exact same angle?
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
    */
  }

  bool Angle::operator<(const Angle &ang) const {
    _sanity_check_(*this);
    _sanity_check_(ang);
    return m_data->unique_id < ang.m_data->unique_id;
  }
  bool Angle::operator>(const Angle &ang) const {
    _sanity_check_(*this);
    _sanity_check_(ang);
    return m_data->unique_id > ang.m_data->unique_id;
  }

  // =======================================================================
  // == CONSTRUCTION =======================================================
  // =======================================================================

  Angle::Impl::Impl(const Atom &a, const Atom &b, const Atom &c,
                    const Molecule &mol)
      : atoms({a, b, c}), molecule(mol) {
  }

  Angle::Angle(const Atom &a, const Atom &b, const Atom &c, const Molecule &mol)
      : m_data(std::make_shared<Impl>(a, b, c, mol)) {
  }

  void Angle::Reset() {
    m_data->atoms.fill(Atom());
    m_data->molecule = Molecule();
    m_data->tag = INT64_MIN;
    m_data->unique_id = INT64_MIN;
    m_data->forcefield_type = FFAngle();
  }

  // =======================================================================
  // == STATE CHECKING =====================================================
  // =======================================================================

  bool Angle::HasType() const {
    _sanity_check_(*this);
    return bool(m_data->forcefield_type);
  }

  // =======================================================================
  // == STATE GETTING ======================================================
  // =======================================================================

  int64_t Angle::GetTag() const {
    _sanity_check_(*this);
    return m_data->tag;
  }

  int64_t Angle::GetID() const {
    _sanity_check_(*this);
    return m_data->unique_id;
  }

  const Molecule &Angle::GetMolecule() const {
    _sanity_check_(*this);
    return m_data->molecule;
  }

  const Angle::AngleAtoms &Angle::GetAtoms() const {
    _sanity_check_(*this);
    return m_data->atoms;
  }

  const FFAngle &Angle::GetType() const {
    _sanity_check_(*this);
    return m_data->forcefield_type;
  }

  int64_t Angle::GetIndex() const {
    _sanity_check_(*this);
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
    _sanity_check_(*this);
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
    _sanity_check_(ang);
    os << "Angle(" << ang.m_data->atoms[0].GetIndex() + 1 << ", "
       << ang.m_data->atoms[1].GetIndex() + 1 << ", "
       << ang.m_data->atoms[2].GetIndex() + 1 << ")";
    return os;
  }
} // namespace indigox
