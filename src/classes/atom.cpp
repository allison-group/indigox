#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/serialise.hpp>

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid atom instance")
#else
#define _sanity_check_(x)
#endif

namespace indigox {

  // =======================================================================
  // == SERIALISATION ======================================================
  // =======================================================================

  template <typename Archive>
  void Atom::Impl::serialise(Archive &archive, const uint32_t) {
    if (INDIGOX_IS_INPUT_ARCHIVE(Archive)) {
      int32_t atomic_number;
      archive(INDIGOX_SERIAL_NVP("element", atomic_number));
      element = GetPeriodicTable()[atomic_number];
    } else {
      archive(INDIGOX_SERIAL_NVP("element", element.GetAtomicNumber()));
    }
    archive(INDIGOX_SERIAL_NVP("molecule", molecule),
            INDIGOX_SERIAL_NVP("formal_charge", formal_charge),
            INDIGOX_SERIAL_NVP("tag", tag),
            INDIGOX_SERIAL_NVP("unique_id", unique_id),
            INDIGOX_SERIAL_NVP("implicit_h_count", implicit_hydrogens),
            INDIGOX_SERIAL_NVP("name", name),
            INDIGOX_SERIAL_NVP("position", position),
            INDIGOX_SERIAL_NVP("partial_charge", partial_charge),
            INDIGOX_SERIAL_NVP("stereochemistry", stereochemistry),
            INDIGOX_SERIAL_NVP("type", forcefield_type),
            INDIGOX_SERIAL_NVP("bonds", bonds),
            INDIGOX_SERIAL_NVP("angles", angles),
            INDIGOX_SERIAL_NVP("dihedrals", dihedrals));
  }

  template <typename Archive>
  void Atom::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(Atom);

  // =======================================================================
  // == OPERATORS ==========================================================
  // =======================================================================

  bool Atom::operator==(const Atom &atm) const {
    return m_data == atm.m_data;
  }

  bool Atom::operator<(const Atom &atm) const {
    _sanity_check_(*this);
    _sanity_check_(atm);
    return m_data->unique_id < atm.m_data->unique_id;
  }

  bool Atom::operator>(const Atom &atm) const {
    _sanity_check_(*this);
    _sanity_check_(atm);
    return m_data->unique_id > atm.m_data->unique_id;
  }

  // =======================================================================
  // == CONSTRUCTION =======================================================
  // =======================================================================

  Atom::Impl::Impl(const Molecule &m, const Element &e, double x, double y,
                   double z, std::string n)
      : molecule(m), element(e), formal_charge(0), charge_group_id(-1),
        residue_id(-1), implicit_hydrogens(0), position(x, y, z), name(n),
        residue_name(""), partial_charge(0.0),
        stereochemistry(AtomStereo::UNDEFINED) {
  }

  Atom::Atom(const Molecule &m, const Element &e, double x, double y, double z,
             std::string n)
      : m_data(std::make_shared<Impl>(m, e, x, y, z, n)) {
  }

  void Atom::Reset() {
    m_data->element = GetPeriodicTable().GetUndefined();
    m_data->molecule = Molecule();
    m_data->formal_charge = INT32_MIN;
    m_data->tag = INT64_MIN;
    m_data->unique_id = INT64_MIN;
    m_data->charge_group_id = -1;
    m_data->implicit_hydrogens = UINT32_MAX;
    m_data->name = "";
    m_data->position = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
    m_data->partial_charge = HUGE_VAL;
    m_data->stereochemistry = AtomStereo::UNDEFINED;
    m_data->forcefield_type = FFAtom();
    m_data->bonds.clear();
    m_data->angles.clear();
    m_data->dihedrals.clear();
  }

  // =======================================================================
  // == STATE CHECKING =====================================================
  // =======================================================================

  int64_t Atom::NumBonds() const {
    _sanity_check_(*this);
    return m_data->bonds.size();
  }

  int64_t Atom::NumAngles() const {
    _sanity_check_(*this);
    return m_data->angles.size();
  }

  int64_t Atom::NumDihedrals() const {
    _sanity_check_(*this);
    return m_data->dihedrals.size();
  }

  bool Atom::HasType() const {
    _sanity_check_(*this);
    return bool(m_data->forcefield_type);
  }

  // =======================================================================
  // == STATE GETTING ======================================================
  // =======================================================================

  const Element &Atom::GetElement() const {
    _sanity_check_(*this);
    return m_data->element;
  }

  int32_t Atom::GetFormalCharge() const {
    _sanity_check_(*this);
    return m_data->formal_charge;
  }

  double Atom::GetPartialCharge() const {
    _sanity_check_(*this);
    return m_data->partial_charge;
  }

  int64_t Atom::GetTag() const {
    _sanity_check_(*this);
    return m_data->tag;
  }

  int64_t Atom::GetID() const {
    _sanity_check_(*this);
    return m_data->unique_id;
  }

  int32_t Atom::GetChargeGroupID() const {
    _sanity_check_(*this);
    return m_data->charge_group_id;
  }

  int32_t Atom::GetResidueID() {
    _sanity_check_(*this);
    m_data->molecule.PerceiveResidues();
    return m_data->residue_id;
  }

  std::string Atom::GetResidueName() {
    _sanity_check_(*this);
    m_data->molecule.PerceiveResidues();
    return m_data->residue_name;
  }

  int32_t Atom::GetImplicitCount() const {
    _sanity_check_(*this);
    return m_data->implicit_hydrogens;
  }

  const Molecule &Atom::GetMolecule() const {
    _sanity_check_(*this);
    return m_data->molecule;
  }

  const std::string &Atom::GetName() const {
    _sanity_check_(*this);
    return m_data->name;
  }

  double Atom::GetX() const {
    _sanity_check_(*this);
    return m_data->position[0];
  }

  double Atom::GetY() const {
    _sanity_check_(*this);
    return m_data->position[1];
  }

  double Atom::GetZ() const {
    _sanity_check_(*this);
    return m_data->position[2];
  }

  AtomStereo Atom::GetStereochemistry() const {
    _sanity_check_(*this);
    return m_data->stereochemistry;
  }

  const Eigen::Vector3d &Atom::GetPosition() const {
    _sanity_check_(*this);
    return m_data->position;
  }

  const Atom::AtomBonds &Atom::GetBonds() const {
    _sanity_check_(*this);
    return m_data->bonds;
  }

  const Atom::AtomAngles &Atom::GetAngles() const {
    _sanity_check_(*this);
    return m_data->angles;
  }

  const Atom::AtomDihedrals &Atom::GetDihedrals() const {
    _sanity_check_(*this);
    return m_data->dihedrals;
  }

  int64_t Atom::GetIndex() const {
    _sanity_check_(*this);
    if (!m_data->molecule) {
      return -1;
    }

    const Molecule::MoleculeAtoms &atoms = m_data->molecule.GetAtoms();
    auto found = std::find(atoms.begin(), atoms.end(), *this);
    return found != atoms.end() ? std::distance(atoms.begin(), found) : -1;
  }

  const FFAtom &Atom::GetType() const {
    _sanity_check_(*this);
    return m_data->forcefield_type;
  }

  // =======================================================================
  // == STATE MODIFYING ====================================================
  // =======================================================================

  int32_t Atom::AddImplicitHydrogen() {
    _sanity_check_(*this);
    return ++m_data->implicit_hydrogens;
  }

  int32_t Atom::RemoveImplicitHydrogen() {
    _sanity_check_(*this);
    return m_data->implicit_hydrogens > 0 ? --m_data->implicit_hydrogens
                                          : m_data->implicit_hydrogens;
  }

  void Atom::AddBond(const Bond &bond) {
    _sanity_check_(*this);
    m_data->bonds.emplace_back(bond);
  }

  void Atom::AddAngle(const Angle &angle) {
    _sanity_check_(*this);
    m_data->angles.emplace_back(angle);
  }

  void Atom::AddDihedral(const Dihedral &dihedral) {
    _sanity_check_(*this);
    m_data->dihedrals.emplace_back(dihedral);
  }

  void Atom::RemoveBond(const Bond &bond) {
    _sanity_check_(*this);
    auto found = std::find(m_data->bonds.begin(), m_data->bonds.end(), bond);
    if (found != m_data->bonds.end()) {
      m_data->bonds.erase(found);
    }
  }

  void Atom::RemoveAngle(const Angle &angle) {
    _sanity_check_(*this);
    auto found = std::find(m_data->angles.begin(), m_data->angles.end(), angle);
    if (found != m_data->angles.end()) {
      m_data->angles.erase(found);
    }
  }

  void Atom::RemoveDihedral(const Dihedral &dihedral) {
    _sanity_check_(*this);
    auto found =
        std::find(m_data->dihedrals.begin(), m_data->dihedrals.end(), dihedral);
    if (found != m_data->dihedrals.end()) {
      m_data->dihedrals.erase(found);
    }
  }

  // =======================================================================
  // == STATE SETTING ======================================================
  // =======================================================================

  void Atom::SetElement(const Element &element) {
    _sanity_check_(*this);
    if (m_data->molecule) {
      m_data->molecule.ModificationMade();
    }
    m_data->element = element;
  }

  void Atom::SetElement(std::string symbol) {
    SetElement(GetPeriodicTable()[symbol]);
  }

  void Atom::SetElement(int32_t atomic_number) {
    SetElement(GetPeriodicTable()[atomic_number]);
  }

  void Atom::SetFormalCharge(int32_t fc) {
    _sanity_check_(*this);
    m_data->formal_charge = fc;
  }

  void Atom::SetPartialCharge(double partial) {
    _sanity_check_(*this);
    m_data->partial_charge = partial;
  }

  void Atom::SetImplicitCount(int32_t count) {
    _sanity_check_(*this);
    count < 0 ? m_data->implicit_hydrogens = 0
              : m_data->implicit_hydrogens = count;
  }

  void Atom::SetTag(int32_t tag) {
    _sanity_check_(*this);
    m_data->tag = tag;
  }

  void Atom::SetChargeGroupID(int32_t id) {
    _sanity_check_(*this);
    m_data->charge_group_id = id;
  }

  void Atom::SetName(std::string name) {
    _sanity_check_(*this);
    m_data->name = name;
  }

  void Atom::SetPosition(double x, double y, double z) {
    _sanity_check_(*this);
    m_data->position = {x, y, z};
  }

  void Atom::SetStereochemistry(AtomStereo stereo) {
    _sanity_check_(*this);
    m_data->stereochemistry = stereo;
  }

  void Atom::SetType(const FFAtom &type) {
    _sanity_check_(*this);
    if (m_data->molecule) {
      if (!m_data->molecule.HasForcefield()) {
        m_data->molecule.SetForcefield(type.GetForcefield());
      } else if (m_data->molecule.GetForcefield() != type.GetForcefield()) {
        throw std::runtime_error(
            "Provided atom type does not match molecule's forcefield");
      }
    }
    m_data->forcefield_type = type;
  }

  // =======================================================================
  // == OUTPUTTING =========================================================
  // =======================================================================

  std::ostream &operator<<(std::ostream &os, const Atom &atm) {
    if (atm) {
      os << "Atom(" << atm.GetIndex() + 1 << ")";
    }
    return os;
  }
} // namespace indigox
