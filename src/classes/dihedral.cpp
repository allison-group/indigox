#include <indigox/classes/atom.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
#include <indigox/utils/serialise.hpp>

#include <memory>
#include <vector>

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid dihedral instance")
#else
#define _sanity_check_(x)
#endif

namespace indigox {

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
  INDIGOX_SERIALISE(Dihedral);

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

  std::ostream &operator<<(std::ostream &os, const Dihedral &dhd) {
    _sanity_check_(dhd);
    os << "Dihedral(" << dhd.m_data->atoms[0].GetIndex() + 1 << ", "
       << dhd.m_data->atoms[1].GetIndex() + 1 << ", " << dhd.m_data->atoms[2].GetIndex() + 1 << ", "
       << dhd.m_data->atoms[3].GetIndex() + 1 << ")";
    return os;
  }
} // namespace indigox
