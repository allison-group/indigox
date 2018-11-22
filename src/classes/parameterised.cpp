#include <boost/math/special_functions/relative_difference.hpp>

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
    wAtom atom;
    TypeCounts types;
    MappedCharge charges;
    MappedAtoms mapped_atoms;
    bool applied;
    
    ParamAtomImpl(Atom& atm) : atom(atm.weak_from_this()), applied(false) { }
  };
  
  // ===========================================================================
  // == ParamAtom Construction/Assignment ======================================
  // ===========================================================================
  
  ParamAtom::ParamAtom() : m_patmdat(nullptr) { }
  ParamAtom::ParamAtom(const ParamAtom& atm) : m_patmdat(atm.m_patmdat) { }
  ParamAtom::ParamAtom(ParamAtom&& atm) : m_patmdat(std::move(atm.m_patmdat)) { }
  ParamAtom& ParamAtom::operator=(const ParamAtom &atm) {
    if (&atm != this) m_patmdat = atm.m_patmdat;
    return *this;
  }
  ParamAtom& ParamAtom::operator=(ParamAtom &&atm) {
    m_patmdat = std::move(atm.m_patmdat);
    return *this;
  }
  ParamAtom::ParamAtom(Atom& atm)
  : m_patmdat(std::make_shared<ParamAtomImpl>(atm)) { }
  
  // ===========================================================================
  // == ParamAtom Data Modification ============================================
  // ===========================================================================
  
  void ParamAtom::MappedWith(Atom& mapped) {
    if (m_patmdat->applied) return;
    if (!mapped.HasType()) throw std::runtime_error("Needs a parameterised atom");
    FFAtom t = mapped.GetType();
    auto t_pos = m_patmdat->types.find(t);
    if (t_pos == m_patmdat->types.end()) m_patmdat->types.emplace(t, 1);
    else ++t_pos->second;
    m_patmdat->charges.push_back(mapped.GetPartialCharge());
    m_patmdat->mapped_atoms.emplace_back(mapped.weak_from_this());
  }
  
  bool ParamAtom::ApplyParameterisation(bool self_consistent) {
    if (m_patmdat->applied) return false;
    if (m_patmdat->charges.empty()) return false;
    double mean = MeanCharge();
    if (self_consistent) {
      if (m_patmdat->types.size() > 1)
        throw std::runtime_error("Types not self-consistent");
      if (std::fabs(mean - MeadianCharge()) > 1e-10)
        throw std::runtime_error("Charges mean/median not equal");
      if (std::fabs(0.0 - StandardDeviationCharge()) > 1e-10)
        throw std::runtime_error("Charge stddev not 0");
    }
    sAtom atm = m_patmdat->atom.lock();
    if (!atm) throw std::runtime_error("Mapped atom missing");
    atm->SetType(GetMostCommonType());
    if (!self_consistent) atm->SetPartialCharge(mean);
    else atm->SetPartialCharge(MeadianCharge());
    m_patmdat->applied = true;
    return true;
  }
  
  // ===========================================================================
  // == ParamAtom Data Retrevial ===============================================
  // ===========================================================================
  
  int64_t ParamAtom::NumSourceAtoms() const {
    return m_patmdat->mapped_atoms.size(); }
  Atom& ParamAtom::GetAtom() const { return *m_patmdat->atom.lock(); }
  double ParamAtom::MeanCharge() const {
    return CalculateMean(m_patmdat->charges.begin(), m_patmdat->charges.end()); }
  double ParamAtom::MeadianCharge() {
    return CalculateMedian(m_patmdat->charges.begin(), m_patmdat->charges.end()); }
  double ParamAtom::StandardDeviationCharge() const {
    return CalculateStandardDeviation(m_patmdat->charges.begin(),
                                      m_patmdat->charges.end()); }
  FFAtom ParamAtom::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    T max = *std::max_element(m_patmdat->types.begin(), m_patmdat->types.end(),
                              [](T& a, T& b) { return a.second < b.second; });
    return max.first; }
  const ParamAtom::TypeCounts& ParamAtom::GetMappedTypeCounts() const {
    return m_patmdat->types; }
  const ParamAtom::MappedCharge& ParamAtom::GetMappedCharges() const {
    return m_patmdat->charges; }
  
  // ===========================================================================
  // == ParamAtom Operators ====================================================
  // ===========================================================================
  
  bool ParamAtom::operator==(const ParamAtom &atm) const {
    return m_patmdat->atom.lock() == atm.m_patmdat->atom.lock(); }
  bool ParamAtom::operator!=(const ParamAtom &atm) const {
    return !(*this == atm); }
  bool ParamAtom::operator<(const ParamAtom &atm) const {
    return (m_patmdat->atom.lock()->GetIndex()
            < atm.m_patmdat->atom.lock()->GetIndex()); }
  bool ParamAtom::operator>(const ParamAtom &atm) const {
    return (m_patmdat->atom.lock()->GetIndex()
            > atm.m_patmdat->atom.lock()->GetIndex()); }
  bool ParamAtom::operator<=(const ParamAtom &atm) const {
    return !(*this > atm); }
  bool ParamAtom::operator>=(const ParamAtom &atm) const {
    return !(*this < atm); }
  ParamAtom::operator bool() const { return bool(m_patmdat); }
  std::ostream& operator<<(std::ostream& os, const ParamAtom& atm) {
    if (atm) os << "Param[" << atm.GetAtom() << "]";
    return os;
  }
  
  // ===========================================================================
  // == ParamBond Data Implementation ==========================================
  // ===========================================================================
  
  struct ParamBond::ParamBondImpl {
    BondAtoms atoms;
    wBond bond;
    TypeCounts types;
    MappedBonds mapped_bonds;
    bool applied;
    
    ParamBondImpl(Atom& a, Atom& b, Bond& bnd)
    : atoms(a.weak_from_this(), b.weak_from_this()), bond(bnd.weak_from_this()),
    applied(false) { }
  };
  
  // ===========================================================================
  // == ParamBond Construction/Assignment ======================================
  // ===========================================================================
  
  ParamBond::ParamBond() : m_pbnddat(nullptr) { }
  ParamBond::ParamBond(const ParamBond& bnd) : m_pbnddat(bnd.m_pbnddat) { }
  ParamBond::ParamBond(ParamBond&& bnd) : m_pbnddat(std::move(bnd.m_pbnddat)) { }
  ParamBond& ParamBond::operator=(const ParamBond &bnd) {
    if (&bnd != this) m_pbnddat = bnd.m_pbnddat;
    return *this;
  }
  ParamBond& ParamBond::operator=(ParamBond &&bnd) {
    m_pbnddat = std::move(bnd.m_pbnddat);
    return *this;
  }
  ParamBond::ParamBond(std::pair<Atom&, Atom&> atms, Bond& bnd)
  : m_pbnddat(std::make_shared<ParamBondImpl>(atms.first, atms.second, bnd)) { }
  
  // ===========================================================================
  // == ParamBond Data Modification ============================================
  // ===========================================================================
  
  void ParamBond::MappedWith(Bond &mapped) {
    if (m_pbnddat->applied) return;
    if (!mapped.HasType()) throw std::runtime_error("Needs a parameterised bond");
    FFBond t = mapped.GetType();
    FFBond t2 = t.GetLinkedType();
    auto t_pos = m_pbnddat->types.find(t);
    auto t2_pos = m_pbnddat->types.find(t2);
    auto end = m_pbnddat->types.end();
    
    if (t_pos == end && t2_pos == end) m_pbnddat->types.emplace(t, 1);
    else if (t2 && t2_pos != end) ++t2_pos->second;
    else ++t_pos->second;
    m_pbnddat->mapped_bonds.emplace_back(mapped.weak_from_this());
  }
  
  bool ParamBond::ApplyParameterisation(bool self_consistent) {
    if (m_pbnddat->applied) return false;
    if (m_pbnddat->mapped_bonds.empty()) return false;
    if (self_consistent && m_pbnddat->types.size() > 1)
      throw std::runtime_error("Types not self-consistent");
    sBond bnd = m_pbnddat->bond.lock();
    if (!bnd) throw std::runtime_error("Mapped bond missing");
    bnd->SetType(GetMostCommonType());
    m_pbnddat->applied = true;
    return true;
  }
  
  // ===========================================================================
  // == ParamBond Data Retrevial ===============================================
  // ===========================================================================
  
  int64_t ParamBond::NumSourceBonds() const {
    return m_pbnddat->mapped_bonds.size(); }
  std::pair<Atom&, Atom&> ParamBond::GetAtoms() const {
    auto atms = m_pbnddat->atoms;
    return {*atms.first.lock(), *atms.second.lock()}; }
  Bond& ParamBond::GetBond() const { return *m_pbnddat->bond.lock(); }
  FFBond ParamBond::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    T max = *std::max_element(m_pbnddat->types.begin(), m_pbnddat->types.end(),
                              [](T& a, T& b) { return a.second < b.second; });
    return max.first; }
  const ParamBond::TypeCounts& ParamBond::GetMappedTypeCounts() const {
    return m_pbnddat->types; }
  
  // ===========================================================================
  // == ParamBond Operators ====================================================
  // ===========================================================================
  
  bool ParamBond::operator==(const ParamBond& bnd) const {
    return m_pbnddat->bond.lock() == bnd.m_pbnddat->bond.lock(); }
  bool ParamBond::operator!=(const ParamBond &bnd) const {
    return !(*this == bnd); }
  bool ParamBond::operator<(const ParamBond &bnd) const {
    return (m_pbnddat->bond.lock()->GetIndex()
            < bnd.m_pbnddat->bond.lock()->GetIndex()); }
  bool ParamBond::operator>(const ParamBond &bnd) const {
    return (m_pbnddat->bond.lock()->GetIndex()
            > bnd.m_pbnddat->bond.lock()->GetIndex()); }
  bool ParamBond::operator<=(const ParamBond &bnd) const {
    return !(*this > bnd); }
  bool ParamBond::operator>=(const ParamBond &bnd) const {
    return !(*this < bnd); }
  ParamBond::operator bool() const { return bool(m_pbnddat); }
  std::ostream& operator<<(std::ostream& os, const ParamBond& bnd) {
    if (bnd) os << "Param[" << bnd.GetBond() << "]";
    return os;
  }
  
  // ===========================================================================
  // == ParamAngle Data Implementation =========================================
  // ===========================================================================
  
  struct ParamAngle::ParamAngleImpl {
    AngleAtoms atoms;
    wAngle angle;
    TypeCounts types;
    MappedAngles mapped_angles;
    bool applied;
    
    ParamAngleImpl(Atom& a, Atom& b, Atom& c, Angle& ang)
    : atoms(a.weak_from_this(), b.weak_from_this(), c.weak_from_this()),
    angle(ang.weak_from_this()), applied(false) { }
  };
  
  // ===========================================================================
  // == ParamAngle Construction/Assignment =====================================
  // ===========================================================================
  
  ParamAngle::ParamAngle() : m_pangdat(nullptr) { }
  ParamAngle::ParamAngle(const ParamAngle& ang) : m_pangdat(ang.m_pangdat) { }
  ParamAngle::ParamAngle(ParamAngle&& ang)
  : m_pangdat(std::move(ang.m_pangdat)) { }
  ParamAngle& ParamAngle::operator=(const ParamAngle &ang) {
    if (&ang != this) m_pangdat = ang.m_pangdat;
    return *this;
  }
  ParamAngle& ParamAngle::operator=(ParamAngle &ang) {
    m_pangdat = std::move(ang.m_pangdat);
    return *this;
  }
  ParamAngle::ParamAngle(stdx::triple<Atom&> atms, Angle& ang)
  : m_pangdat(std::make_shared<ParamAngleImpl>(atms.first, atms.second,
                                               atms.third, ang)) { }
  
  // ===========================================================================
  // == ParamAngle Data Modification ===========================================
  // ===========================================================================
  
  void ParamAngle::MappedWith(Angle& mapped) {
    if (m_pangdat->applied) return;
    if (!mapped.HasType()) throw std::runtime_error("Needs a parameterised angle");
    FFAngle t = mapped.GetType();
    FFAngle t2 = t.GetLinkedType();
    auto t_pos = m_pangdat->types.find(t);
    auto t2_pos = m_pangdat->types.find(t2);
    auto end = m_pangdat->types.end();
    
    if (t_pos == end && t2_pos == end) m_pangdat->types.emplace(t, 1);
    else if (t2 && t2_pos != end) ++t2_pos->second;
    else ++t_pos->second;
    
    m_pangdat->mapped_angles.emplace_back(mapped.weak_from_this());
  }
  
  bool ParamAngle::ApplyParameterisation(bool self_consistent) {
    if (m_pangdat->mapped_angles.empty()) return false;
    if (m_pangdat->applied) return false;
    if (self_consistent && m_pangdat->types.size() > 1)
      throw std::runtime_error("Types not self-consistent");
    sAngle ang = m_pangdat->angle.lock();
    if (!ang) throw std::runtime_error("Mapped angle missing");
    ang->SetType(GetMostCommonType());
    m_pangdat->applied = true;
    return true;
  }
  
  // ===========================================================================
  // == ParamAngle Data Retrevial ==============================================
  // ===========================================================================
  
  int64_t ParamAngle::NumSourceAngles() const {
    return m_pangdat->mapped_angles.size(); }
  stdx::triple<Atom&> ParamAngle::GetAtoms() const {
    auto atms = m_pangdat->atoms;
    return {*atms.first.lock(), *atms.second.lock(), *atms.third.lock()}; }
  Angle& ParamAngle::GetAngle() const { return *m_pangdat->angle.lock(); }
  FFAngle ParamAngle::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    T max = *std::max_element(m_pangdat->types.begin(), m_pangdat->types.end(),
                              [](T& a, T& b) { return a.second < b.second; });
    return max.first; }
  const ParamAngle::TypeCounts& ParamAngle::GetMappedTypeCounts() const {
    return m_pangdat->types; }
  
  // ===========================================================================
  // == ParamAngle Operators ===================================================
  // ===========================================================================
  
  bool ParamAngle::operator==(const ParamAngle& ang) const {
    return m_pangdat->angle.lock() == ang.m_pangdat->angle.lock(); }
  bool ParamAngle::operator!=(const ParamAngle &ang) const {
    return !(*this == ang); }
  bool ParamAngle::operator<(const ParamAngle &ang) const {
    return (m_pangdat->angle.lock()->GetIndex()
            < ang.m_pangdat->angle.lock()->GetIndex()); }
  bool ParamAngle::operator>(const ParamAngle &ang) const {
    return (m_pangdat->angle.lock()->GetIndex()
            > ang.m_pangdat->angle.lock()->GetIndex()); }
  bool ParamAngle::operator<=(const ParamAngle &ang) const {
    return !(*this > ang); }
  bool ParamAngle::operator>=(const ParamAngle &ang) const {
    return !(*this < ang); }
  ParamAngle::operator bool() const { return bool(m_pangdat); }
  std::ostream& operator<<(std::ostream& os, const ParamAngle& ang) {
    if (ang) os << "Param[" << ang.GetAngle() << "]";
    return os;
  }
  
  // ===========================================================================
  // == ParamDihedral Data Implementation ======================================
  // ===========================================================================
  
  struct ParamDihedral::ParamDihedralImpl {
    DihedralAtoms atoms;
    wDihedral dihedral;
    TypeCounts types;
    MappedDihedrals mapped_dihedrals;
    bool applied;
    
    ParamDihedralImpl(Atom& a, Atom& b, Atom& c, Atom& d, Dihedral& dhd)
    : atoms(a.weak_from_this(), b.weak_from_this(), c.weak_from_this(),
            d.weak_from_this()), dihedral(dhd.weak_from_this()), applied(false)
    { }
  };
  
  // ===========================================================================
  // == ParamDihedral Construction/Assignment ==================================
  // ===========================================================================
  
  ParamDihedral::ParamDihedral() : m_pdhddat(nullptr) { }
  ParamDihedral::ParamDihedral(const ParamDihedral& dhd)
  : m_pdhddat(dhd.m_pdhddat) { }
  ParamDihedral::ParamDihedral(ParamDihedral&& dhd)
  : m_pdhddat(std::move(dhd.m_pdhddat)) { }
  ParamDihedral& ParamDihedral::operator=(const ParamDihedral &dhd) {
    if (&dhd != this) m_pdhddat = dhd.m_pdhddat;
    return *this;
  }
  ParamDihedral& ParamDihedral::operator=(ParamDihedral &&dhd) {
    m_pdhddat = std::move(dhd.m_pdhddat);
    return *this;
  }
  ParamDihedral::ParamDihedral(stdx::quad<Atom&> atms, Dihedral& dhd)
  : m_pdhddat(std::make_shared<ParamDihedralImpl>(atms.first, atms.second,
                                                  atms.third, atms.fourth, dhd))
  { }
  
  // ===========================================================================
  // == ParamDihedral Data Modification ========================================
  // ===========================================================================
  
  void ParamDihedral::MappedWith(Dihedral& mapped) {
    if (m_pdhddat->applied) return;
    // Can't throw with dihedrals as they can have no types assigned
    if (!mapped.HasType()) return;
    
    TypeGroup t = mapped.GetTypes();
    auto t_pos = m_pdhddat->types.emplace(t, 1);
    if (!t_pos.second) ++(t_pos.first->second);
    m_pdhddat->mapped_dihedrals.emplace_back(mapped.weak_from_this());
  }
  
  bool ParamDihedral::ApplyParameterisation(bool self_consistent) {
    if (m_pdhddat->mapped_dihedrals.empty()) return false;
    if (m_pdhddat->applied) return false;
    if (self_consistent && m_pdhddat->types.size() > 1)
      throw std::runtime_error("Types not self-consistent");
    sDihedral dhd = m_pdhddat->dihedral.lock();
    if (!dhd) throw std::runtime_error("Mapped dihedral missing");
    dhd->SetTypes(GetMostCommonType());
    m_pdhddat->applied = true;
    return true;
  }
  
  // ===========================================================================
  // == ParamDihedral Data Retrevial ===========================================
  // ===========================================================================
  
  int64_t ParamDihedral::NumSourceDihedral() const {
    return m_pdhddat->mapped_dihedrals.size(); }
  stdx::quad<Atom&> ParamDihedral::GetParameterisedAtoms() const {
    return {*m_pdhddat->atoms.first.lock(), *m_pdhddat->atoms.second.lock(),
      *m_pdhddat->atoms.third.lock(), *m_pdhddat->atoms.fourth.lock()};
  }
  Dihedral& ParamDihedral::GetDihedral() const {
    return *m_pdhddat->dihedral.lock(); }
  ParamDihedral::TypeGroup ParamDihedral::GetMostCommonType() const {
    using T = TypeCounts::value_type;
    T max = *std::max_element(m_pdhddat->types.begin(), m_pdhddat->types.end(),
                              [](T& a, T& b) { return a.second < b.second; });
    return max.first; }
  const ParamDihedral::TypeCounts& ParamDihedral::GetMappedTypeCounts() const {
    return m_pdhddat->types;
  }
  
  // ===========================================================================
  // == ParamDihedral Operators ================================================
  // ===========================================================================
  
  bool ParamDihedral::operator==(const ParamDihedral& ang) const {
    return m_pdhddat->dihedral.lock() == ang.m_pdhddat->dihedral.lock(); }
  bool ParamDihedral::operator!=(const ParamDihedral &ang) const {
    return !(*this == ang); }
  bool ParamDihedral::operator<(const ParamDihedral &ang) const {
    return (m_pdhddat->dihedral.lock()->GetIndex()
            < ang.m_pdhddat->dihedral.lock()->GetIndex()); }
  bool ParamDihedral::operator>(const ParamDihedral &ang) const {
    return (m_pdhddat->dihedral.lock()->GetIndex()
            > ang.m_pdhddat->dihedral.lock()->GetIndex()); }
  bool ParamDihedral::operator<=(const ParamDihedral &ang) const {
    return !(*this > ang); }
  bool ParamDihedral::operator>=(const ParamDihedral &ang) const {
    return !(*this < ang); }
  ParamDihedral::operator bool() const { return bool(m_pdhddat); }
  std::ostream& operator<<(std::ostream& os, const ParamDihedral& ang) {
    if (ang) os << "Param[" << ang.GetDihedral() << "]";
    return os;
  }
  
  // ===========================================================================
  // == ParamMolecule Data Implementation ======================================
  // ===========================================================================
  struct ParamMolecule::ParamMoleculeImpl {
    sMolecule mol;
    ParamAtoms atoms;
    ParamBonds bonds;
    ParamAngles angles;
    ParamDihedrals dihedrals;
    std::vector<ParamAtom> nonsc_atoms;
    
    ParamMoleculeImpl(Molecule& m) : mol(m.shared_from_this()) {
      for (sAtom atm : m.GetAtoms()) atoms.emplace(atm, ParamAtom(*atm));
      for (sBond bnd : m.GetBonds()) {
        auto as = bnd->GetAtoms();
        sAtom a = as.first.shared_from_this();
        sAtom b = as.second.shared_from_this();
        if (a > b) a.swap(b);
        bonds.emplace(std::make_pair(a, b), ParamBond(as, *bnd));
      }
      for (sAngle ang : m.GetAngles()) {
        auto as = ang->GetAtoms();
        sAtom a = as.first.shared_from_this();
        sAtom b = as.second.shared_from_this();
        sAtom c = as.third.shared_from_this();
        if (a > c) a.swap(c);
        angles.emplace(stdx::make_triple(a, b, c), ParamAngle(as, *ang));
      }
      for (sDihedral dhd : m.GetDihedrals()) {
        auto as = dhd->GetAtoms();
        sAtom a = as.first.shared_from_this();
        sAtom b = as.second.shared_from_this();
        sAtom c = as.third.shared_from_this();
        sAtom d = as.fourth.shared_from_this();
        if (a > d) {
          a.swap(d);
          b.swap(c);
        }
        dihedrals.emplace(stdx::make_quad(a, b, c, d), ParamDihedral(as, *dhd));
      }
    }
  };
  
  // ===========================================================================
  // == ParamMolecule Construction/Assignment ==================================
  // ===========================================================================
  
  ParamMolecule::ParamMolecule() : m_pmoldat(nullptr) { }
  ParamMolecule::ParamMolecule(const ParamMolecule& pmol)
  : m_pmoldat(pmol.m_pmoldat) { }
  ParamMolecule::ParamMolecule(ParamMolecule&& pmol)
  : m_pmoldat(std::move(pmol.m_pmoldat)) { }
  ParamMolecule& ParamMolecule::operator=(const ParamMolecule& pmol) {
    if (&pmol != this) m_pmoldat = pmol.m_pmoldat;
    return *this;
  }
  ParamMolecule& ParamMolecule::operator=(ParamMolecule&& pmol) {
    m_pmoldat = std::move(pmol.m_pmoldat);
    return *this;
  }
  ParamMolecule::ParamMolecule(Molecule& mol)
  : m_pmoldat(std::make_shared<ParamMoleculeImpl>(mol)) { }
  
  // ===========================================================================
  // == ParamMolecule Data Modification ========================================
  // ===========================================================================
  
  void ParamMolecule::ApplyParameteristion(bool sc) {
    for (auto& atm : m_pmoldat->atoms) {
      bool param = atm.second.ApplyParameterisation(sc);
      if (!sc && param) m_pmoldat->nonsc_atoms.emplace_back(atm.second);
    }
    for (auto& bnd : m_pmoldat->bonds) bnd.second.ApplyParameterisation(sc);
    for (auto& ang : m_pmoldat->angles) ang.second.ApplyParameterisation(sc);
    for (auto& dhd : m_pmoldat->dihedrals) dhd.second.ApplyParameterisation(sc);
  }
  
  // ===========================================================================
  // == ParamMolecule Data Retrevial ===========================================
  // ===========================================================================
  
  ParamAtom ParamMolecule::GetAtom(Atom &atm) const {
    auto pos = m_pmoldat->atoms.find(atm.shared_from_this());
    return pos == m_pmoldat->atoms.end() ? ParamAtom() : pos->second;
  }
  
  ParamBond ParamMolecule::GetBond(Bond& bnd) const {
    sAtom a = bnd.GetSourceAtom().shared_from_this();
    sAtom b = bnd.GetTargetAtom().shared_from_this();
    return GetBond(std::make_pair(a, b));
  }
  
  ParamBond ParamMolecule::GetBond(PBond atms) const {
    if (atms.first > atms.second) atms.first.swap(atms.second);
    auto pos = m_pmoldat->bonds.find(atms);
    return pos == m_pmoldat->bonds.end() ? ParamBond() : pos->second;
  }
  
  ParamAngle ParamMolecule::GetAngle(Angle &ang) const {
    auto atms = ang.GetAtoms();
    sAtom a = atms.first.shared_from_this();
    sAtom b = atms.second.shared_from_this();
    sAtom c = atms.third.shared_from_this();
    return GetAngle(stdx::make_triple(a, b, c));
  }
  
  ParamAngle ParamMolecule::GetAngle(PAngle atms) const {
    if (atms.first > atms.third) atms.first.swap(atms.third);
    auto pos = m_pmoldat->angles.find(atms);
    return pos == m_pmoldat->angles.end() ? ParamAngle() : pos->second;
  }
  
  ParamDihedral ParamMolecule::GetDihedral(Dihedral &dhd) {
    auto atms = dhd.GetAtoms();
    sAtom a = atms.first.shared_from_this();
    sAtom b = atms.second.shared_from_this();
    sAtom c = atms.third.shared_from_this();
    sAtom d = atms.fourth.shared_from_this();
    return GetDihedral(stdx::make_quad(a, b, c, d));
  }
  
  ParamDihedral ParamMolecule::GetDihedral(PDihedral atms) {
    if (atms.fourth < atms.first) {
      atms.first.swap(atms.fourth);
      atms.second.swap(atms.third);
    }
    auto pos = m_pmoldat->dihedrals.find(atms);
    if (pos != m_pmoldat->dihedrals.end()) return pos->second;
    Dihedral& newD = m_pmoldat->mol->NewDihedral(*atms.first, *atms.second,
                                                *atms.third, *atms.fourth);
    ParamDihedral newPD(newD.GetAtoms(), newD);
    m_pmoldat->dihedrals.emplace(atms, newPD);
    return newPD;
  }
  
  std::vector<ParamAtom> ParamMolecule::GetAtoms() const {
    std::vector<ParamAtom> atms; atms.reserve(m_pmoldat->atoms.size());
    for (auto& a : m_pmoldat->atoms) atms.emplace_back(a.second);
    std::sort(atms.begin(), atms.end(),
              [](ParamAtom& a, ParamAtom& b) {
                return a.GetAtom().GetIndex() < b.GetAtom().GetIndex();
              });
    return atms;
  }
  
  std::vector<ParamBond> ParamMolecule::GetBonds() const {
    std::vector<ParamBond> atms; atms.reserve(m_pmoldat->bonds.size());
    for (auto& a : m_pmoldat->bonds) atms.emplace_back(a.second);
    std::sort(atms.begin(), atms.end(),
              [](ParamBond& a, ParamBond& b) {
                size_t aa = a.GetAtoms().first.GetIndex();
                size_t ab = a.GetAtoms().second.GetIndex();
                if (aa > ab) std::swap(aa, ab);
                size_t ba = b.GetAtoms().first.GetIndex();
                size_t bb = b.GetAtoms().second.GetIndex();
                if (ba > bb) std::swap(ba, bb);
                return (aa == ba) ? (ab < bb) : (aa < ba);
              });
    return atms;
  }
  
  std::vector<ParamAngle> ParamMolecule::GetAngles() const {
    std::vector<ParamAngle> atms; atms.reserve(m_pmoldat->angles.size());
    for (auto& a : m_pmoldat->angles) atms.emplace_back(a.second);
    std::sort(atms.begin(), atms.end(),
              [](ParamAngle& a, ParamAngle& b) {
                size_t aa = a.GetAtoms().first.GetIndex();
                size_t ab = a.GetAtoms().second.GetIndex();
                size_t ac = a.GetAtoms().third.GetIndex();
                if (aa > ac) std::swap(aa, ac);
                size_t ba = b.GetAtoms().first.GetIndex();
                size_t bb = b.GetAtoms().second.GetIndex();
                size_t bc = b.GetAtoms().third.GetIndex();
                if (ba > bc) std::swap(ba, bc);
                if (ab == bb && aa == ba) return ac < bc;
                else if (ab == bb) return aa < ba;
                else return ab < bb;
              });
    return atms;
  }
  
  std::vector<ParamDihedral> ParamMolecule::GetDihedrals() const {
    std::vector<ParamDihedral> atms; atms.reserve(m_pmoldat->dihedrals.size());
    for (auto& a : m_pmoldat->dihedrals) atms.emplace_back(a.second);
    std::sort(atms.begin(), atms.end(),
              [](ParamDihedral& a, ParamDihedral& b) {
                size_t aa = a.GetParameterisedAtoms().first.GetIndex();
                size_t ab = a.GetParameterisedAtoms().second.GetIndex();
                size_t ac = a.GetParameterisedAtoms().third.GetIndex();
                size_t ad = a.GetParameterisedAtoms().fourth.GetIndex();
                if (ab > ac) { std::swap(ab, ac); std::swap(aa, ad); }
                size_t ba = b.GetParameterisedAtoms().first.GetIndex();
                size_t bb = b.GetParameterisedAtoms().second.GetIndex();
                size_t bc = b.GetParameterisedAtoms().third.GetIndex();
                size_t bd = b.GetParameterisedAtoms().fourth.GetIndex();
                if (bb > bc) { std::swap(bb, bc); std::swap(ba, bd); }
                if (ab == bb && ac == bc && aa == ba) return ad < bd;
                else if (ab == bb && ac == bc) return aa < ba;
                else if (ab == bb) return ac < bc;
                else return ab < bb;
              });
    return atms;
  }
  
  // ===========================================================================
  // == ParamMolecule Operators ================================================
  // ===========================================================================
  
}
