#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>

namespace indigox {

  // ===========================================================================
  // == ParamAtom Data Implementation ==========================================
  // ===========================================================================

  struct ParamAtom::ParamAtomImpl {
    Atom atom;
    TypeCounts types;
    MappedCharge charges;
    bool applied;

    ParamAtomImpl(const Atom &atm) : atom(atm), applied(false) {
    }
  };

  // ===========================================================================
  // == ParamAtom Construction/Assignment ======================================
  // ===========================================================================

  ParamAtom::ParamAtom(const Atom &atm)
      : m_data(std::make_shared<ParamAtomImpl>(atm)) {
  }

  // ===========================================================================
  // == ParamAtom Data Modification ============================================
  // ===========================================================================

  void ParamAtom::MappedWith(const Atom &mapped) {
    if (m_data->applied)
      return;
    if (!mapped.HasType())
      throw std::runtime_error("Needs a parameterised atom");
    FFAtom t = mapped.GetType();
    auto t_pos = m_data->types.find(t);
    if (t_pos == m_data->types.end())
      m_data->types.emplace(t, 1);
    else
      ++t_pos->second;
    m_data->charges.push_back(mapped.GetPartialCharge());
  }

  bool ParamAtom::ApplyParameterisation(bool self_consistent) {
    if (m_data->applied)
      return false;
    if (m_data->charges.empty())
      return false;
    double mean = MeanCharge();
    if (self_consistent) {
      if (m_data->types.size() > 1)
        throw std::runtime_error("Types not self-consistent");
      if (std::fabs(mean - MeadianCharge()) > 1e-10)
        throw std::runtime_error("Charges mean/median not equal");
      if (std::fabs(0.0 - StandardDeviationCharge()) > 1e-10)
        throw std::runtime_error("Charge stddev not 0");
    }
    if (!m_data->atom)
      throw std::runtime_error("Mapped atom missing");
    m_data->atom.SetType(GetMostCommonType());
    if (!self_consistent)
      m_data->atom.SetPartialCharge(mean);
    else
      m_data->atom.SetPartialCharge(MeadianCharge());
    m_data->applied = true;
    return true;
  }

  // ===========================================================================
  // == ParamAtom Data Retrevial ===============================================
  // ===========================================================================

  int64_t ParamAtom::NumSourceAtoms() const {
    return m_data->charges.size();
  }
  const Atom &ParamAtom::GetAtom() const {
    return m_data->atom;
  }
  double ParamAtom::MeanCharge() const {
    return CalculateMean(m_data->charges.begin(), m_data->charges.end());
  }
  double ParamAtom::MeadianCharge() {
    return CalculateMedian(m_data->charges.begin(), m_data->charges.end());
  }
  double ParamAtom::StandardDeviationCharge() const {
    return CalculateStandardDeviation(m_data->charges.begin(),
                                      m_data->charges.end());
  }
  const FFAtom &ParamAtom::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    return std::max_element(m_data->types.begin(), m_data->types.end(),
                            [](T &a, T &b) { return a.second < b.second; })
        ->first;
  }
  const ParamAtom::TypeCounts &ParamAtom::GetMappedTypeCounts() const {
    return m_data->types;
  }
  const ParamAtom::MappedCharge &ParamAtom::GetMappedCharges() const {
    return m_data->charges;
  }

  // ===========================================================================
  // == ParamAtom Operators ====================================================
  // ===========================================================================

  bool ParamAtom::operator==(const ParamAtom &atm) const {
    return m_data->atom == atm.m_data->atom;
  }
  bool ParamAtom::operator<(const ParamAtom &atm) const {
    return (m_data->atom.GetIndex() < atm.m_data->atom.GetIndex());
  }
  bool ParamAtom::operator>(const ParamAtom &atm) const {
    return (m_data->atom.GetIndex() > atm.m_data->atom.GetIndex());
  }
  std::ostream &operator<<(std::ostream &os, const ParamAtom &atm) {
    if (atm)
      os << "Param[" << atm.GetAtom() << "]";
    return os;
  }

  // ===========================================================================
  // == ParamBond Data Implementation ==========================================
  // ===========================================================================

  struct ParamBond::ParamBondImpl {
    BondAtoms atoms;
    Bond bond;
    TypeCounts types;
    bool applied;

    ParamBondImpl(BondAtoms atms, const Bond &bnd)
        : atoms(atms), bond(bnd), applied(false) {
    }
  };

  // ===========================================================================
  // == ParamBond Construction/Assignment ======================================
  // ===========================================================================

  ParamBond::ParamBond(BondAtoms atms, const Bond &bnd)
      : m_data(std::make_shared<ParamBondImpl>(atms, bnd)) {
  }

  // ===========================================================================
  // == ParamBond Data Modification ============================================
  // ===========================================================================

  void ParamBond::MappedWith(const Bond &mapped) {
    if (m_data->applied)
      return;
    if (!mapped.HasType())
      throw std::runtime_error("Needs a parameterised bond");
    FFBond t = mapped.GetType();
    FFBond t2 = t.GetLinkedType();
    auto t_pos = m_data->types.find(t);
    auto t2_pos = m_data->types.find(t2);
    auto end = m_data->types.end();

    if (t_pos == end && t2_pos == end)
      m_data->types.emplace(t, 1);
    else if (t2 && t2_pos != end)
      ++t2_pos->second;
    else
      ++t_pos->second;
  }

  bool ParamBond::ApplyParameterisation(bool self_consistent) {
    if (m_data->applied)
      return false;
    if (m_data->types.empty())
      return false;
    if (self_consistent && m_data->types.size() > 1)
      throw std::runtime_error("Types not self-consistent");
    if (!m_data->bond)
      throw std::runtime_error("Mapped bond missing");
    m_data->bond.SetType(GetMostCommonType());
    m_data->applied = true;
    return true;
  }

  // ===========================================================================
  // == ParamBond Data Retrevial ===============================================
  // ===========================================================================

  int64_t ParamBond::NumSourceBonds() const {
    int64_t count = 0;
    for (auto &countable : m_data->types)
      count += countable.second;
    return count;
  }

  const ParamBond::BondAtoms &ParamBond::GetAtoms() const {
    return m_data->atoms;
  }

  const Bond &ParamBond::GetBond() const {
    return m_data->bond;
  }

  const FFBond &ParamBond::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    return std::max_element(m_data->types.begin(), m_data->types.end(),
                            [](T &a, T &b) { return a.second < b.second; })
        ->first;
  }

  const ParamBond::TypeCounts &ParamBond::GetMappedTypeCounts() const {
    return m_data->types;
  }

  // ===========================================================================
  // == ParamBond Operators ====================================================
  // ===========================================================================

  bool ParamBond::operator==(const ParamBond &bnd) const {
    return m_data->bond == bnd.m_data->bond;
  }

  bool ParamBond::operator<(const ParamBond &bnd) const {
    return (m_data->bond.GetIndex() < bnd.m_data->bond.GetIndex());
  }
  bool ParamBond::operator>(const ParamBond &bnd) const {
    return (m_data->bond.GetIndex() > bnd.m_data->bond.GetIndex());
  }
  std::ostream &operator<<(std::ostream &os, const ParamBond &bnd) {
    if (bnd)
      os << "Param[" << bnd.GetBond() << "]";
    return os;
  }

  // ===========================================================================
  // == ParamAngle Data Implementation =========================================
  // ===========================================================================

  struct ParamAngle::ParamAngleImpl {
    AngleAtoms atoms;
    Angle angle;
    TypeCounts types;
    bool applied;

    ParamAngleImpl(AngleAtoms atms, const Angle &ang)
        : atoms(atms), angle(ang), applied(false) {
    }
  };

  // ===========================================================================
  // == ParamAngle Construction/Assignment =====================================
  // ===========================================================================

  ParamAngle::ParamAngle(AngleAtoms atms, const Angle &ang)
      : m_data(std::make_shared<ParamAngleImpl>(atms, ang)) {
  }

  // ===========================================================================
  // == ParamAngle Data Modification ===========================================
  // ===========================================================================

  void ParamAngle::MappedWith(const Angle &mapped) {
    if (m_data->applied)
      return;
    if (!mapped.HasType())
      throw std::runtime_error("Needs a parameterised angle");
    FFAngle t = mapped.GetType();
    FFAngle t2 = t.GetLinkedType();
    auto t_pos = m_data->types.find(t);
    auto t2_pos = m_data->types.find(t2);
    auto end = m_data->types.end();

    if (t_pos == end && t2_pos == end)
      m_data->types.emplace(t, 1);
    else if (t2 && t2_pos != end)
      ++t2_pos->second;
    else
      ++t_pos->second;
  }

  bool ParamAngle::ApplyParameterisation(bool self_consistent) {
    if (m_data->types.empty())
      return false;
    if (m_data->applied)
      return false;
    if (self_consistent && m_data->types.size() > 1)
      throw std::runtime_error("Types not self-consistent");

    if (!m_data->angle)
      throw std::runtime_error("Mapped angle missing");
    m_data->angle.SetType(GetMostCommonType());
    m_data->applied = true;
    return true;
  }

  // ===========================================================================
  // == ParamAngle Data Retrevial ==============================================
  // ===========================================================================

  int64_t ParamAngle::NumSourceAngles() const {
    int64_t count = 0;
    for (auto &countable : m_data->types)
      count += countable.second;
    return count;
  }

  const ParamAngle::AngleAtoms &ParamAngle::GetAtoms() const {
    return m_data->atoms;
  }

  const Angle &ParamAngle::GetAngle() const {
    return m_data->angle;
  }

  const FFAngle &ParamAngle::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    return std::max_element(m_data->types.begin(), m_data->types.end(),
                            [](T &a, T &b) { return a.second < b.second; })
        ->first;
  }

  const ParamAngle::TypeCounts &ParamAngle::GetMappedTypeCounts() const {
    return m_data->types;
  }

  // ===========================================================================
  // == ParamAngle Operators ===================================================
  // ===========================================================================

  bool ParamAngle::operator==(const ParamAngle &ang) const {
    return m_data->angle == ang.m_data->angle;
  }

  bool ParamAngle::operator<(const ParamAngle &ang) const {
    return (m_data->angle.GetIndex() < ang.m_data->angle.GetIndex());
  }
  bool ParamAngle::operator>(const ParamAngle &ang) const {
    return (m_data->angle.GetIndex() > ang.m_data->angle.GetIndex());
  }
  std::ostream &operator<<(std::ostream &os, const ParamAngle &ang) {
    if (ang)
      os << "Param[" << ang.GetAngle() << "]";
    return os;
  }

  // ===========================================================================
  // == ParamDihedral Data Implementation ======================================
  // ===========================================================================

  struct ParamDihedral::ParamDihedralImpl {
    DihedralAtoms atoms;
    Dihedral dihedral;
    TypeCounts types;
    bool applied;

    ParamDihedralImpl(DihedralAtoms atms, const Dihedral &dhd)
        : atoms(atms), dihedral(dhd), applied(false) {
    }
  };

  // ===========================================================================
  // == ParamDihedral Construction/Assignment ==================================
  // ===========================================================================

  ParamDihedral::ParamDihedral(DihedralAtoms atms, const Dihedral &dhd)
      : m_data(std::make_shared<ParamDihedralImpl>(atms, dhd)) {
  }

  // ===========================================================================
  // == ParamDihedral Data Modification ========================================
  // ===========================================================================

  void ParamDihedral::MappedWith(const Dihedral &mapped) {
    if (m_data->applied)
      return;
    // Can't throw with dihedrals as they can have no types assigned
    if (!mapped.HasType())
      return;
    // Only map with same priority dihedrals.
    //  To avoid GROMOS duplication in eg NH3 group.
    if (mapped.GetPriority() >= 0 &&
        mapped.GetPriority() != m_data->dihedral.GetPriority())
      return;

    TypeGroup t = mapped.GetTypes();
    auto t_pos = m_data->types.emplace(t, 1);
    if (!t_pos.second)
      ++(t_pos.first->second);
  }

  bool ParamDihedral::ApplyParameterisation(bool self_consistent) {
    if (m_data->types.empty())
      return false;
    if (m_data->applied)
      return false;
    if (self_consistent && m_data->types.size() > 1)
      throw std::runtime_error("Types not self-consistent");

    if (!m_data->dihedral)
      throw std::runtime_error("Mapped dihedral missing");
    m_data->dihedral.SetTypes(GetMostCommonType());
    m_data->applied = true;
    return true;
  }

  // ===========================================================================
  // == ParamDihedral Data Retrevial ===========================================
  // ===========================================================================

  int64_t ParamDihedral::NumSourceDihedral() const {
    int64_t count = 0;
    for (auto &countable : m_data->types)
      count += countable.second;
    return count;
  }

  const ParamDihedral::DihedralAtoms &ParamDihedral::GetAtoms() const {
    return m_data->atoms;
  }
  const Dihedral &ParamDihedral::GetDihedral() const {
    return m_data->dihedral;
  }
  const ParamDihedral::TypeGroup &ParamDihedral::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    return std::max_element(m_data->types.begin(), m_data->types.end(),
                            [](T &a, T &b) { return a.second < b.second; })
        ->first;
  }

  const ParamDihedral::TypeCounts &ParamDihedral::GetMappedTypeCounts() const {
    return m_data->types;
  }

  // ===========================================================================
  // == ParamDihedral Operators ================================================
  // ===========================================================================

  bool ParamDihedral::operator==(const ParamDihedral &ang) const {
    return m_data->dihedral == ang.m_data->dihedral;
  }

  bool ParamDihedral::operator<(const ParamDihedral &ang) const {
    return (m_data->dihedral.GetIndex() < ang.m_data->dihedral.GetIndex());
  }

  bool ParamDihedral::operator>(const ParamDihedral &ang) const {
    return (m_data->dihedral.GetIndex() > ang.m_data->dihedral.GetIndex());
  }
  std::ostream &operator<<(std::ostream &os, const ParamDihedral &ang) {
    if (ang)
      os << "Param[" << ang.GetDihedral() << "]";
    return os;
  }

  // ===========================================================================
  // == ParamMolecule Data Implementation ======================================
  // ===========================================================================
  struct ParamMolecule::ParamMoleculeImpl {
    Molecule mol;
    std::vector<ParamAtom> atoms;
    std::vector<ParamBond> bonds;
    std::vector<ParamAngle> angles;
    std::vector<ParamDihedral> dihedrals;

    ParamAtoms atom_indices;
    ParamBonds bond_indices;
    ParamAngles angle_indices;
    ParamDihedrals dihedral_indices;
    std::vector<ParamAtom> nonsc_atoms;

    ParamMoleculeImpl(const Molecule &m) : mol(m) {
      for (const Atom &atm : mol.GetAtoms()) {
        atoms.emplace_back(atm);
        atom_indices.emplace(atm, atom_indices.size());
      }

      for (const Bond &bnd : mol.GetBonds()) {
        Atom a = bnd.GetAtoms()[0];
        Atom b = bnd.GetAtoms()[1];
        if (a > b) {
          std::swap(a, b);
        }
        bonds.emplace_back(std::make_pair(a, b), bnd);
        bond_indices.emplace(std::make_pair(a, b), bond_indices.size());
      }

      for (const Angle &ang : mol.GetAngles()) {
        Atom a = ang.GetAtoms()[0];
        Atom b = ang.GetAtoms()[1];
        Atom c = ang.GetAtoms()[2];
        if (a > c) {
          std::swap(a, c);
        }
        angles.emplace_back(stdx::make_triple(a, b, c), ang);
        angle_indices.emplace(stdx::make_triple(a, b, c), angle_indices.size());
      }

      for (const Dihedral &dhd : mol.GetDihedrals()) {
        Atom a = dhd.GetAtoms()[0];
        Atom b = dhd.GetAtoms()[1];
        Atom c = dhd.GetAtoms()[2];
        Atom d = dhd.GetAtoms()[3];
        if (a > d) {
          std::swap(a, d);
          std::swap(b, c);
        }
        dihedrals.emplace_back(stdx::make_quad(a, b, c, d), dhd);
        dihedral_indices.emplace(stdx::make_quad(a, b, c, d),
                                 dihedral_indices.size());
      }
    }
  };

  // ===========================================================================
  // == ParamMolecule Construction/Assignment ==================================
  // ===========================================================================

  ParamMolecule::ParamMolecule(const Molecule &mol)
      : m_data(std::make_shared<ParamMoleculeImpl>(mol)) {
  }

  // ===========================================================================
  // == ParamMolecule Data Modification ========================================
  // ===========================================================================

  void ParamMolecule::ApplyParameteristion(bool sc) {
    for (ParamAtom atm : m_data->atoms) {
      bool param = atm.ApplyParameterisation(sc);
      if (!sc && param)
        m_data->nonsc_atoms.emplace_back(atm);
    }
    for (ParamBond bnd : m_data->bonds)
      bnd.ApplyParameterisation(sc);
    for (ParamAngle ang : m_data->angles)
      ang.ApplyParameterisation(sc);
    for (ParamDihedral dhd : m_data->dihedrals)
      dhd.ApplyParameterisation(sc);
  }

  // ===========================================================================
  // == ParamMolecule Data Retrevial ===========================================
  // ===========================================================================

  const ParamAtom &ParamMolecule::GetAtom(const Atom &atm) const {
    return m_data->atoms[m_data->atom_indices.find(atm)->second];
  }

  const ParamBond &ParamMolecule::GetBond(const Bond &bnd) const {
    return GetBond(bnd.GetAtoms()[0], bnd.GetAtoms()[1]);
  }

  const ParamBond &ParamMolecule::GetBond(const Atom &a, const Atom &b) const {
    PBond dat = (a > b) ? std::make_pair(b, a) : std::make_pair(a, b);
    return m_data->bonds[m_data->bond_indices.find(dat)->second];
  }

  const ParamAngle &ParamMolecule::GetAngle(const Angle &ang) const {
    return GetAngle(ang.GetAtoms()[0], ang.GetAtoms()[1], ang.GetAtoms()[2]);
  }

  const ParamAngle &ParamMolecule::GetAngle(const Atom &a, const Atom &b,
                                            const Atom &c) const {
    PAngle dat =
        (a > c) ? stdx::make_triple(c, b, a) : stdx::make_triple(a, b, c);
    return m_data->angles[m_data->angle_indices.find(dat)->second];
  }

  const ParamDihedral &ParamMolecule::GetDihedral(const Dihedral &dhd) {
    return GetDihedral(dhd.GetAtoms()[0], dhd.GetAtoms()[1], dhd.GetAtoms()[2],
                       dhd.GetAtoms()[3]);
  }

  const ParamDihedral &ParamMolecule::GetDihedral(const Atom &a, const Atom &b,
                                                  const Atom &c,
                                                  const Atom &d) {
    PDihedral dat =
        (a > d) ? stdx::make_quad(d, c, b, a) : stdx::make_quad(a, b, c, d);
    auto pos = m_data->dihedral_indices.find(dat);
    if (pos != m_data->dihedral_indices.end())
      return m_data->dihedrals[pos->second];

    Dihedral newD = m_data->mol.NewDihedral(a, b, c, d);
    m_data->dihedrals.emplace_back(dat, newD);
    m_data->dihedral_indices.emplace(dat, m_data->dihedral_indices.size());
    return m_data->dihedrals.back();
  }

  const std::vector<ParamAtom> &ParamMolecule::GetAtoms() const {
    return m_data->atoms;
  }

  const std::vector<ParamBond> &ParamMolecule::GetBonds() const {
    return m_data->bonds;
  }

  const std::vector<ParamAngle> &ParamMolecule::GetAngles() const {
    return m_data->angles;
  }

  const std::vector<ParamDihedral> &ParamMolecule::GetDihedrals() const {
    return m_data->dihedrals;
  }

  // ===========================================================================
  // == ParamMolecule Operators ================================================
  // ===========================================================================

} // namespace indigox
