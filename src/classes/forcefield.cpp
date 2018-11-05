//#include <algorithm>
#include <array>
#include <initializer_list>

#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/quad.hpp>
#include <indigox/utils/triple.hpp>

#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/forcefield_test.hpp>

namespace indigox {

  // ===========================================================================
  // == FFAtom Data Implementation =============================================
  // ===========================================================================
  
  struct FFAtom::FFAtomImpl {
    int32_t id;
    std::string name;
    int32_t element;
    uint32_t implicit_hydrogens;
    Forcefield ff;
    
    FFAtomImpl() = default;
    FFAtomImpl(int32_t i, std::string n, const Element& e, uint32_t h,
               const Forcefield& f)
    : id(i), name(n), element(e.GetAtomicNumber()), implicit_hydrogens(h), ff(f)
    { }
    
    template <class Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("id", id),
              INDIGOX_SERIAL_NVP("name", name),
              INDIGOX_SERIAL_NVP("element", element),
              INDIGOX_SERIAL_NVP("implicitH", implicit_hydrogens),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };
  
  // ===========================================================================
  // == FFAtom Construction and Assignment =====================================
  // ===========================================================================
  
  FFAtom::FFAtom() : m_ffatmdat(nullptr) { }
  FFAtom::FFAtom(const FFAtom& atm) : m_ffatmdat(atm.m_ffatmdat) { }
  FFAtom::FFAtom(FFAtom&& atm) : m_ffatmdat(std::move(atm.m_ffatmdat)) { }
  FFAtom& FFAtom::operator=(const FFAtom &atm) {
    if (&atm != this) m_ffatmdat = atm.m_ffatmdat;
    return *this;
  }
  FFAtom& FFAtom::operator=(FFAtom &&atm) {
    m_ffatmdat = std::move(atm.m_ffatmdat);
    return *this;
  }
  FFAtom::FFAtom(int32_t id, std::string name, const Element& e, uint32_t h,
                 const Forcefield& ff)
  : m_ffatmdat(std::make_shared<FFAtomImpl>(id, name, e, h, ff)) { }
  
  // ===========================================================================
  // == FFAtom Serialisation ===================================================
  // ===========================================================================
  
  template <class Archive>
  void FFAtom::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_ffatmdat));
  }
  INDIGOX_SERIALISE(FFAtom);
  
  // ===========================================================================
  // == FFAtom Data Retrevial ==================================================
  // ===========================================================================
  
  int32_t FFAtom::GetID() const { return m_ffatmdat->id; }
  std::string FFAtom::GetName() const { return m_ffatmdat->name; }
  Forcefield& FFAtom::GetForcefield() const { return m_ffatmdat->ff; }
  uint32_t FFAtom::GetImplicitHydrogenCount() const {
    return m_ffatmdat->implicit_hydrogens; }
  Element FFAtom::GetElement() const {
    return GetPeriodicTable()[m_ffatmdat->element]; }
  
  // ===========================================================================
  // == FFAtom Operators =======================================================
  // ===========================================================================
  
  bool FFAtom::operator==(const FFAtom &atm) const {
    return GetID() == atm.GetID() && GetName() == atm.GetName();
  }
  bool FFAtom::operator!=(const FFAtom &atm) const { return !(*this == atm); }
  bool FFAtom::operator<(const FFAtom &atm) const {
    return GetID() < atm.GetID();
  }
  bool FFAtom::operator>(const FFAtom &atm) const {
    return GetID() > atm.GetID();
  }
  bool FFAtom::operator<=(const FFAtom &atm) const { return !(*this > atm); }
  bool FFAtom::operator>=(const FFAtom &atm) const { return !(*this < atm); }
  FFAtom::operator bool() const { return bool(m_ffatmdat); }
  
  std::ostream& operator<<(std::ostream& os, const FFAtom& atm) {
    if (atm) os << "FFAtom(" << atm.GetID() << ")";
    return os;
  }
  
  // ===========================================================================
  // == FFBond Data implementation =============================================
  // ===========================================================================
  
  struct FFBond::FFBondImpl {
    BondType type;
    int32_t id;
    FFBond::DataStore raw_data;
    FFBond::AllowedMask allowed_parameters;
    FFBond linked_bond;
    Forcefield ff;
    
    FFBondImpl() = default;
    FFBondImpl(BondType t, int32_t i, FFParam data, const Forcefield& f)
    : type(t), id(i), linked_bond(), ff(f) {
      if (t == BondType::Harmonic) allowed_parameters.from_uint32(3);
      if (t == BondType::Quartic) allowed_parameters.from_uint32(3);
      raw_data.fill(0);
      if (data.size() != allowed_parameters.count())
        throw std::range_error("Incorrect item count");
      std::copy(data.begin(), data.end(), raw_data.begin());
    }
    
    template <class Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("type", type),
              INDIGOX_SERIAL_NVP("id", id),
              INDIGOX_SERIAL_NVP("raw_data", raw_data),
              INDIGOX_SERIAL_NVP("allowed", allowed_parameters),
              INDIGOX_SERIAL_NVP("linked_bond_type", linked_bond),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };
  
  // ===========================================================================
  // == FFBond Construction and Assignment =====================================
  // ===========================================================================
  
  FFBond::FFBond() : m_ffbnddat(nullptr) { }
  FFBond::FFBond(const FFBond& bnd) : m_ffbnddat(bnd.m_ffbnddat) { }
  FFBond::FFBond(FFBond&& bnd) : m_ffbnddat(std::move(bnd.m_ffbnddat)) { }
  FFBond& FFBond::operator=(const FFBond &bnd) {
    if (&bnd != this) m_ffbnddat = bnd.m_ffbnddat;
    return *this;
  }
  FFBond& FFBond::operator=(FFBond &&bnd) {
    m_ffbnddat = std::move(bnd.m_ffbnddat);
    return *this;
  }
  FFBond::FFBond(BondType type, int32_t id, FFParam params, const Forcefield& ff)
  : m_ffbnddat(std::make_shared<FFBondImpl>(type, id, params, ff)) { }
  
  // ===========================================================================
  // == FFBond Serialisation ===================================================
  // ===========================================================================
  
  template <class Archive>
  void FFBond::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_ffbnddat));
  }
  INDIGOX_SERIALISE(FFBond)
  
  // ===========================================================================
  // == FFBond Data Retrevial ==================================================
  // ===========================================================================
  
  double FFBond::GetForceConstant() const {
    return (m_ffbnddat->allowed_parameters[Allow_ForceConstant]
            ? m_ffbnddat->raw_data[Store_ForceConstant]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFBond::GetIdealLength() const {
    return (m_ffbnddat->allowed_parameters[Allow_IdealLength]
            ? m_ffbnddat->raw_data[Store_IdealLength]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  BondType FFBond::GetType() const { return m_ffbnddat->type; }
  int32_t FFBond::GetID() const { return m_ffbnddat->id; }
  FFBond& FFBond::GetLinkedType() const { return m_ffbnddat->linked_bond; }
  Forcefield& FFBond::GetForcefield() const { return m_ffbnddat->ff; }
  
  // ===========================================================================
  // == FFBond Operators =======================================================
  // ===========================================================================
  
  bool FFBond::operator==(const FFBond& bnd) const {
    return (GetType() == bnd.GetType() && GetID() == bnd.GetID()
            && m_ffbnddat->raw_data == bnd.m_ffbnddat->raw_data);
  }
  bool FFBond::operator!=(const FFBond &bnd) const { return !(*this == bnd); }
  bool FFBond::operator<(const FFBond &bnd) const {
    if (GetType() < bnd.GetType()) return true;
    if (GetID() < bnd.GetID()) return true;
    return false;
  }
  bool FFBond::operator>(const FFBond &bnd) const {
    if (GetType() > bnd.GetType()) return true;
    if (GetID() > bnd.GetID()) return true;
    return false;
  }
  bool FFBond::operator<=(const FFBond &bnd) const { return !(*this > bnd); }
  bool FFBond::operator>=(const FFBond &bnd) const { return !(*this < bnd); }
  FFBond::operator bool() const { return bool(m_ffbnddat); }
  
  std::ostream& operator<<(std::ostream& os, BondType type) {
    switch (type) {
      case FFBond::Type::Empty: return os << "Empty";
      case FFBond::Type::Harmonic: return os << "Harmonic";
      case FFBond::Type::Quartic: return os << "Quartic";
      default: return os;
    }
  }
  
  std::ostream& operator<<(std::ostream& os, const FFBond& bnd) {
    if (bnd) os << "FFBond(" << bnd.GetType() << ", " << bnd.GetID() << ")";
    return os;
  }
  
  // ===========================================================================
  // == FFAngle Data Implementation ============================================
  // ===========================================================================
  
  struct FFAngle::FFAngleImpl {
    AngleType type;
    int32_t id;
    FFAngle::DataStore raw_data;
    FFAngle::AllowedMask allowed_parameters;
    FFAngle linked_angle;
    Forcefield ff;
    
    FFAngleImpl() = default;
    FFAngleImpl(AngleType t, int32_t i, FFParam data, const Forcefield& f)
    : type(t), id(i), linked_angle(), ff(f) {
      if (t == AngleType::Harmonic) allowed_parameters.from_uint32(3);
      if (t == AngleType::CosineHarmonic) allowed_parameters.from_uint32(3);
      raw_data.fill(0);
      if (data.size() != allowed_parameters.count())
        throw std::range_error("Incorrect item count");
      std::copy(data.begin(), data.end(), raw_data.begin());
    }
    
    template <class Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("type", type),
              INDIGOX_SERIAL_NVP("id", id),
              INDIGOX_SERIAL_NVP("raw_data", raw_data),
              INDIGOX_SERIAL_NVP("allwed_mask", allowed_parameters),
              INDIGOX_SERIAL_NVP("link", linked_angle),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };
  
  // ===========================================================================
  // == FFAngle Construction and Assignment ====================================
  // ===========================================================================
  
  FFAngle::FFAngle() : m_ffangdat(nullptr) { }
  FFAngle::FFAngle(const FFAngle& ang) : m_ffangdat(ang.m_ffangdat) { }
  FFAngle::FFAngle(FFAngle&& ang) : m_ffangdat(std::move(ang.m_ffangdat)) { }
  FFAngle& FFAngle::operator=(const FFAngle& ang) {
    if (&ang != this) m_ffangdat = ang.m_ffangdat;
    return *this;
  }
  FFAngle& FFAngle::operator=(FFAngle &&ang) {
    m_ffangdat = std::move(ang.m_ffangdat);
    return *this;
  }
  FFAngle::FFAngle(AngleType type, int32_t id, FFParam par, const Forcefield& ff)
  : m_ffangdat(std::make_shared<FFAngleImpl>(type, id, par, ff)) { }
  
  // ===========================================================================
  // == FFAngle Serialisation ==================================================
  // ===========================================================================
  
  template <class Archive>
  void FFAngle::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_ffangdat));
  }
  INDIGOX_SERIALISE(FFAngle);
  
  // ===========================================================================
  // == FFAngle Data Retrevial =================================================
  // ===========================================================================
  
  double FFAngle::GetForceConstant() const {
    return (m_ffangdat->allowed_parameters[Allow_ForceConstant]
            ? m_ffangdat->raw_data[Store_ForceConstant]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFAngle::GetIdealAngle() const {
    return (m_ffangdat->allowed_parameters[Allow_IdealAngle]
            ? m_ffangdat->raw_data[Store_IdealAngle]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  AngleType FFAngle::GetType() const { return m_ffangdat->type; }
  int32_t FFAngle::GetID() const { return m_ffangdat->id; }
  FFAngle& FFAngle::GetLinkedType() const { return m_ffangdat->linked_angle; }
  Forcefield& FFAngle::GetForcefield() const { return m_ffangdat->ff; }
  
  // ===========================================================================
  // == FFAngle Operators ======================================================
  // ===========================================================================
  
  bool FFAngle::operator==(const FFAngle& bnd) const {
    return (GetType() == bnd.GetType() && GetID() == bnd.GetID()
            && m_ffangdat->raw_data == bnd.m_ffangdat->raw_data);
  }
  bool FFAngle::operator!=(const FFAngle &bnd) const { return !(*this == bnd); }
  bool FFAngle::operator<(const FFAngle &bnd) const {
    if (GetType() < bnd.GetType()) return true;
    if (GetID() < bnd.GetID()) return true;
    return false;
  }
  bool FFAngle::operator>(const FFAngle &bnd) const {
    if (GetType() > bnd.GetType()) return true;
    if (GetID() > bnd.GetID()) return true;
    return false;
  }
  bool FFAngle::operator<=(const FFAngle &bnd) const { return !(*this > bnd); }
  bool FFAngle::operator>=(const FFAngle &bnd) const { return !(*this < bnd); }
  FFAngle::operator bool() const { return bool(m_ffangdat); }
  
  std::ostream& operator<<(std::ostream& os, AngleType type) {
    switch (type) {
      case FFAngle::Type::Empty: return os << "Empty";
      case FFAngle::Type::Harmonic: return os << "Harmonic";
      case FFAngle::Type::CosineHarmonic: return os << "CosineHarmonic";
      default: return os;
    }
  }
  
  std::ostream& operator<<(std::ostream& os, const FFAngle& ang) {
    if (ang) os << "FFAngle(" << ang.GetType() << ", " << ang.GetID() << ")";
    return os;
  }
  
  // ===========================================================================
  // == FFDihedral Data Implementation =========================================
  // ===========================================================================
  
  struct FFDihedral::FFDihedralImpl {
    DihedralType type;
    int32_t id;
    FFDihedral::DataStore raw_data;
    FFDihedral::AllowedMask allowed_parameters;
    Forcefield ff;
    
    FFDihedralImpl() = default;
    FFDihedralImpl(DihedralType t, int32_t i, FFParam data, const Forcefield& f)
    : type(t), id(i), ff(f) {
      if (t == DihedralType::Proper) allowed_parameters.from_uint32(7);
      if (t == DihedralType::Improper) allowed_parameters.from_uint32(10);
      raw_data.fill(0);
      if (data.size() != allowed_parameters.count())
        throw std::range_error("Incorrect item count");
      std::copy(data.begin(), data.end(), raw_data.begin());
    }
    
    template <class Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("type", type),
              INDIGOX_SERIAL_NVP("id", id),
              INDIGOX_SERIAL_NVP("raw_data", raw_data),
              INDIGOX_SERIAL_NVP("allwed_mask", allowed_parameters),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };
  
  // ===========================================================================
  // == FFDihedral Construction and Assignment ====================================
  // ===========================================================================
  
  FFDihedral::FFDihedral() : m_ffdhddat(nullptr) { }
  FFDihedral::FFDihedral(const FFDihedral& ang) : m_ffdhddat(ang.m_ffdhddat) { }
  FFDihedral::FFDihedral(FFDihedral&& ang)
  : m_ffdhddat(std::move(ang.m_ffdhddat)) { }
  FFDihedral& FFDihedral::operator=(const FFDihedral& ang) {
    if (&ang != this) m_ffdhddat = ang.m_ffdhddat;
    return *this;
  }
  FFDihedral& FFDihedral::operator=(FFDihedral &&ang) {
    m_ffdhddat = std::move(ang.m_ffdhddat);
    return *this;
  }
  FFDihedral::FFDihedral(DihedralType type, int32_t id, FFParam par,
                         const Forcefield& ff)
  : m_ffdhddat(std::make_shared<FFDihedralImpl>(type, id, par, ff)) { }
  
  // ===========================================================================
  // == FFDihedral Serialisation ===============================================
  // ===========================================================================
  
  template <class Archive>
  void FFDihedral::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_ffdhddat));
  }
  INDIGOX_SERIALISE(FFDihedral);
  
  // ===========================================================================
  // == FFDihedral Data Retrevial ==============================================
  // ===========================================================================
  
  double FFDihedral::GetPhaseShift() const {
    return (m_ffdhddat->allowed_parameters[Allow_PhaseShift]
            ? m_ffdhddat->raw_data[Store_PhaseShift]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFDihedral::GetForceConstant() const {
    return (m_ffdhddat->allowed_parameters[Allow_ForceConstant]
            ? m_ffdhddat->raw_data[Store_ForceConstant]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  int32_t FFDihedral::GetMultiplicity() const {
    return (m_ffdhddat->allowed_parameters[Allow_Multiplicity]
            ? static_cast<int32_t>(m_ffdhddat->raw_data[Store_Multiplicity])
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFDihedral::GetIdealAngle() const {
    return (m_ffdhddat->allowed_parameters[Allow_IdealAngle]
            ? m_ffdhddat->raw_data[Store_IdealAngle]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  DihedralType FFDihedral::GetType() const { return m_ffdhddat->type; }
  int32_t FFDihedral::GetID() const { return m_ffdhddat->id; }
  Forcefield& FFDihedral::GetForcefield() const { return m_ffdhddat->ff; }
  
  // ===========================================================================
  // == FFDihedral Operators ===================================================
  // ===========================================================================
  
  bool FFDihedral::operator==(const FFDihedral& bnd) const {
    return (GetType() == bnd.GetType() && GetID() == bnd.GetID()
            && m_ffdhddat->raw_data == bnd.m_ffdhddat->raw_data);
  }
  bool FFDihedral::operator!=(const FFDihedral &bnd) const {
    return !(*this == bnd);
  }
  bool FFDihedral::operator<(const FFDihedral &bnd) const {
    if (GetType() < bnd.GetType()) return true;
    if (GetID() < bnd.GetID()) return true;
    return false;
  }
  bool FFDihedral::operator>(const FFDihedral &bnd) const {
    if (GetType() > bnd.GetType()) return true;
    if (GetID() > bnd.GetID()) return true;
    return false;
  }
  bool FFDihedral::operator<=(const FFDihedral &bnd) const {
    return !(*this > bnd);
  }
  bool FFDihedral::operator>=(const FFDihedral &bnd) const {
    return !(*this < bnd);
  }
  FFDihedral::operator bool() const { return bool(m_ffdhddat); }
  
  std::ostream& operator<<(std::ostream& os, DihedralType dhd) {
    switch (dhd) {
      case FFDihedral::Type::Empty: return (os << "Empty");
      case FFDihedral::Type::Proper: return (os << "Proper");
      case FFDihedral::Type::Improper: return (os << "Improper");
      default: return os;
    }
  }
  
  std::ostream& operator<<(std::ostream& os, const FFDihedral& dhd) {
    if (dhd) os << "FFDihedral(" << dhd.GetType() << ", " << dhd.GetID() << ")";
    return os;
  }
  
  // ===========================================================================
  // == Forcefield Data Implementation =========================================
  // ===========================================================================
  
  struct Forcefield::ForcefieldImpl {
    FFFamily family;
    std::string name;
    Forcefield::AtomTypes atoms;
    Forcefield::BondTypes bonds;
    Forcefield::AngleTypes angles;
    Forcefield::DihedralTypes dihedrals;
    
    ForcefieldImpl() = default;
    ForcefieldImpl(FFFamily fam, std::string nme) : family(fam), name(nme) {
      using BndStr = Forcefield::BondTypes::mapped_type;
      using AngStr = Forcefield::AngleTypes::mapped_type;
      using DhdStr = Forcefield::DihedralTypes::mapped_type;
      
      if (family == FFFamily::GROMOS) {
        bonds.emplace(BondType::Harmonic, BndStr());
        bonds.emplace(BondType::Quartic, BndStr());
        angles.emplace(AngleType::Harmonic, AngStr());
        angles.emplace(AngleType::CosineHarmonic, AngStr());
        dihedrals.emplace(DihedralType::Proper, DhdStr());
        dihedrals.emplace(DihedralType::Improper, DhdStr());
      }
      if (family == FFFamily::Other) {
        bonds.emplace(BondType::Harmonic, BndStr());
        bonds.emplace(BondType::Quartic, BndStr());
        bonds.emplace(BondType::Morse, BndStr());
        bonds.emplace(BondType::Cubic, BndStr());
        bonds.emplace(BondType::FENE, BndStr());
        angles.emplace(AngleType::Harmonic, AngStr());
        angles.emplace(AngleType::CosineHarmonic, AngStr());
        angles.emplace(AngleType::UreyBradley, AngStr());
        angles.emplace(AngleType::Quartic, AngStr());
        dihedrals.emplace(DihedralType::Proper, DhdStr());
        dihedrals.emplace(DihedralType::Improper, DhdStr());
        dihedrals.emplace(DihedralType::RyckaertBellemans, DhdStr());
        dihedrals.emplace(DihedralType::PeriodicImproper, DhdStr());
        dihedrals.emplace(DihedralType::Fourier, DhdStr());
        dihedrals.emplace(DihedralType::Restricted, DhdStr());
      }
    }
    
    template <class Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("family", family),
              INDIGOX_SERIAL_NVP("name", name),
              INDIGOX_SERIAL_NVP("atoms", atoms),
              INDIGOX_SERIAL_NVP("bonds", bonds),
              INDIGOX_SERIAL_NVP("angles", angles),
              INDIGOX_SERIAL_NVP("dihedrals", dihedrals));
    }
    
    bool Contains(BondType t) { return bonds.find(t) != bonds.end(); }
    bool Contains(BondType t, int32_t i) {
      if (!Contains(t)) return false;
      return (std::find_if(bonds.at(t).begin(), bonds.at(t).end(),
                          [&i](auto& b) { return b.GetID() == i;})
              != bonds.at(t).end());
    }
    bool Contains(AngleType t) { return angles.find(t) != angles.end(); }
    bool Contains(AngleType t, int32_t i) {
      if (!Contains(t)) return false;
      return (std::find_if(angles.at(t).begin(), angles.at(t).end(),
                           [&i](auto& b) { return b.GetID() == i;})
              != angles.at(t).end());
    }
    bool Contains(DihedralType t) { return dihedrals.find(t) != dihedrals.end(); }
    bool Contains(DihedralType t, int32_t i) {
      if (!Contains(t)) return false;
      return (std::find_if(dihedrals.at(t).begin(), dihedrals.at(t).end(),
                           [&i](auto& b) { return b.GetID() == i;})
              != dihedrals.at(t).end());
    }
  };
  
  // ===========================================================================
  // == Forcefield Construction and Assignment =================================
  // ===========================================================================
  
  Forcefield::Forcefield() : m_ffdat(nullptr) { }
  Forcefield::Forcefield(const Forcefield& ff) : m_ffdat(ff.m_ffdat) { }
  Forcefield::Forcefield(Forcefield&& ff) : m_ffdat(std::move(ff.m_ffdat)) { }
  Forcefield& Forcefield::operator=(const Forcefield& ff) {
    if (&ff != this) m_ffdat = ff.m_ffdat;
    return *this;
  }
  Forcefield& Forcefield::operator=(Forcefield &&ff) {
    m_ffdat = std::move(ff.m_ffdat);
    return *this;
  }
  Forcefield::Forcefield(FFFamily family, std::string name)
  : m_ffdat(std::make_shared<ForcefieldImpl>(family, name)) { }
  
  // ===========================================================================
  // == Forcefield Serialisation ===============================================
  // ===========================================================================
  
  template <class Archive>
  void Forcefield::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_ffdat));
  }
  INDIGOX_SERIALISE(Forcefield);
  
  // ===========================================================================
  // == Forcefield Data Modification ===========================================
  // ===========================================================================
  
  FFBond& Forcefield::NewBondType(BondType type, int32_t id, FFParam param) {
    if (!m_ffdat->Contains(type))
      throw std::out_of_range("Unsupported bond type");
    if (m_ffdat->Contains(type, id))
      throw std::out_of_range("Bond ID already exists");
    m_ffdat->bonds[type].emplace_back(FFBond(type, id, param, *this));
    return m_ffdat->bonds[type].back();
  }
  
  FFAngle& Forcefield::NewAngleType(AngleType type, int32_t id, FFParam param) {
    if (!m_ffdat->Contains(type))
      throw std::out_of_range("Unsupported angle type");
    if (m_ffdat->Contains(type, id))
      throw std::out_of_range("Angle ID already exists");
    m_ffdat->angles[type].emplace_back(FFAngle(type, id, param, *this));
    return m_ffdat->angles[type].back();
  }
  
  FFDihedral& Forcefield::NewDihedralType(DihedralType type, int32_t id,
                                          FFParam param) {
    if (!m_ffdat->Contains(type))
      throw std::out_of_range("Unsupported dihedral type");
    if (m_ffdat->Contains(type, id))
      throw std::out_of_range("Dihedral ID already exists");
    m_ffdat->dihedrals[type].emplace_back(FFDihedral(type, id, param, *this));
    return m_ffdat->dihedrals[type].back();
  }
  
  FFAtom& Forcefield::NewAtomType(int id, std::string name, const Element& element,
                                  uint32_t implictH) {
    FFAtom atm(id, name, element, implictH, *this);
    if (std::find(m_ffdat->atoms.begin(), m_ffdat->atoms.end(), atm)
        != m_ffdat->atoms.end())
      throw std::out_of_range("Atom type already exists");
    m_ffdat->atoms.emplace_back(atm);
    return m_ffdat->atoms.back();
  }
  
  void Forcefield::ReserveAtomTypes(size_t sz) { m_ffdat->atoms.reserve(sz); }
  
  FFBond& Forcefield::NewBondType(BondType t, int32_t i, double a, double b) {
    return NewBondType(t, i, {a, b});
  }
  
  void Forcefield::LinkBondTypes(FFBond &a, FFBond &b) {
    a.m_ffbnddat->linked_bond = b;
    b.m_ffbnddat->linked_bond = a;
  }
  
  void Forcefield::ReserveBondTypes(BondType type, size_t size) {
    if (m_ffdat->Contains(type)) m_ffdat->bonds.at(type).reserve(size);
  }
  
  FFAngle& Forcefield::NewAngleType(AngleType t, int32_t i, double a, double b) {
    return NewAngleType(t, i, {a, b});
  }
  
  void Forcefield::LinkAngleTypes(FFAngle &a, FFAngle &b) {
    a.m_ffangdat->linked_angle = b;
    b.m_ffangdat->linked_angle = a;
  }
  
  void Forcefield::ReserveAngleTypes(AngleType type, size_t size) {
    if (m_ffdat->Contains(type)) m_ffdat->angles.at(type).reserve(size);
  }
  
  FFDihedral& Forcefield::NewDihedralType(DihedralType type, int32_t id,
                                          double a, double b, double c) {
    return NewDihedralType(type, id, {b,a,c});
  }
  
  FFDihedral& Forcefield::NewDihedralType(DihedralType type, int32_t id,
                                          double a, double b) {
    return NewDihedralType(type, id, {b,a});
  }
  
  void Forcefield::ReserveDihedralTypes(DihedralType type, size_t size) {
    if (m_ffdat->Contains(type)) m_ffdat->dihedrals.at(type).reserve(size);
  }
  
// ============================================================================
// == Forcefield Data Retrevial ===============================================
// ============================================================================
  
  template <typename InputIt>
  size_t sum_of_sizes(InputIt begin, InputIt end) {
    return std::accumulate(begin, end, 0, [](size_t a, auto& b) {
      return a + b.second.size(); });
  }
  FFAtom& Forcefield::GetAtomType(std::string name) const {
    auto bgn = m_ffdat->atoms.begin();
    auto end = m_ffdat->atoms.end();
    auto fnd = [&name](FFAtom& atm) { return atm.GetName() == name; };
    auto pos = std::find_if(bgn, end, fnd);
    if (pos == end) throw std::out_of_range("Atom name does not exist");
    return *pos;
  }
  FFAtom& Forcefield::GetAtomType(int32_t id) const {
    auto bgn = m_ffdat->atoms.begin();
    auto end = m_ffdat->atoms.end();
    auto fnd = [&id](auto& Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    if (pos == end) throw std::out_of_range("Atom ID does not exist");
    return *pos;
  }
  size_t Forcefield::NumAtomTypes() const { return m_ffdat->atoms.size(); }
  FFBond& Forcefield::GetBondType(BondType type, int32_t id) const {
    if (!m_ffdat->Contains(type))
      throw std::out_of_range("Unsupported bond type");
    auto bgn = m_ffdat->bonds.at(type).begin();
    auto end = m_ffdat->bonds.at(type).end();
    auto fnd = [&id](auto& Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    if (pos == end) throw std::out_of_range("Bond ID does not exist");
    return *pos;
  }
  FFBond& Forcefield::GetBondType(int32_t id) const {
    for (auto& types: m_ffdat->bonds) {
      try { return GetBondType(types.first, id); }
      catch (const std::out_of_range&) { continue; }
    }
    throw std::out_of_range("Bond ID does not exist");
  }
  size_t Forcefield::NumBondTypes() const {
    return sum_of_sizes(m_ffdat->bonds.begin(), m_ffdat->bonds.end());
  }
  size_t Forcefield::NumBondTypes(BondType type) const {
    return m_ffdat->Contains(type) ? m_ffdat->bonds.at(type).size() : 0;
  }
  FFAngle& Forcefield::GetAngleType(AngleType type, int32_t id) const {
    if (!m_ffdat->Contains(type))
      throw std::out_of_range("Unsupported angle type");
    auto bgn = m_ffdat->angles.at(type).begin();
    auto end = m_ffdat->angles.at(type).end();
    auto fnd = [&id](auto& Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    if (pos == end) throw std::out_of_range("Angle ID does not exist");
    return *pos;
  }
  FFAngle& Forcefield::GetAngleType(int id) const {
    for (auto& types: m_ffdat->angles) {
      try { return GetAngleType(types.first, id); }
      catch (const std::out_of_range&) { continue; }
    }
    throw std::out_of_range("Angle ID does not exist");
  }
  size_t Forcefield::NumAngleTypes() const {
    return sum_of_sizes(m_ffdat->angles.begin(), m_ffdat->angles.end());
  }
  size_t Forcefield::NumAngleTypes(AngleType type) const {
    return m_ffdat->Contains(type) ? m_ffdat->angles.at(type).size() : 0;
  }
  FFDihedral& Forcefield::GetDihedralType(DihedralType type, int32_t id) const {
    if (!m_ffdat->Contains(type))
      throw std::out_of_range("Unsupported dihedral type");
    auto bgn = m_ffdat->dihedrals.at(type).begin();
    auto end = m_ffdat->dihedrals.at(type).end();
    auto fnd = [&id](auto& Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    if (pos == end) throw std::out_of_range("Dihedral ID does not exist");
    return *pos;
  }
  FFDihedral& Forcefield::GetDihedralType(int32_t id) const {
    for (auto& types: m_ffdat->dihedrals) {
      try { return GetDihedralType(types.first, id); }
      catch (const std::out_of_range&) { continue; }
    }
    throw std::out_of_range("Dihedral ID does not exist");
  }
  size_t Forcefield::NumDihedralTypes() const {
    return sum_of_sizes(m_ffdat->dihedrals.begin(), m_ffdat->dihedrals.end());
  }
  size_t Forcefield::NumDihedralTypes(DihedralType type) const {
    return m_ffdat->Contains(type) ? m_ffdat->dihedrals.at(type).size() : 0;
  }
  FFFamily Forcefield::GetFamily() const { return m_ffdat->family; }
  std::string Forcefield::GetName() const { return m_ffdat->name; }
  
  // ===========================================================================
  // == Forcefield Operators ===================================================
  // ===========================================================================
  
  bool Forcefield::operator==(const Forcefield &ff) const {
    return m_ffdat == ff.m_ffdat;
  }
  bool Forcefield::operator!=(const Forcefield &ff) const {
    return !(*this == ff);
  }
  Forcefield::operator bool() const { return bool(m_ffdat); }
  
  std::ostream& operator<<(std::ostream& os, const Forcefield& ff) {
    if (ff) os << "Forcefield(" << ff.GetName() << ")";
    return os;
  }
  
  // ===========================================================================
  // == Hardcoded forcefield generations =======================================
  // ===========================================================================
  Forcefield GenerateGROMOS54A7() {
    Forcefield ff(FFFamily::GROMOS, "54A7");
    // Add atom types
    PeriodicTable PT = GetPeriodicTable();
    std::vector<stdx::quad<int, std::string, int, int>> atom_dat = {
      {2, "OM", 8, 0}, {1, "O", 8, 0}, {4, "OE", 8, 0}, {39, "CChl", 6, 0},
      {32, "F", 9, 0}, {52, "OUrea", 8, 0}, {15, "CH2", 6, 2}, {20, "HC", 1, 0},
      {37, "NA+", 11, 0}, {33, "CL", 17, 0}, {44, "ODmso", 8, 0},
      {40, "CLChl", 17, 0}, {45, "CCl4", 6, 0}, {31, "AR", 18, 0},
      {38, "CL-", 17, 0}, {35, "CMet", 6, 0}, {21, "H", 1, 0},
      {27, "ZN2+", 30, 0}, {28, "MG2+", 12, 0}, {46, "CLCl4", 17, 0},
      {36, "OMet", 8, 0}, {54, "CH3p", 6, 3}, {51, "CUrea", 6, 0},
      {30, "P,SI", 15, 0}, {7, "NT", 7, 0}, {26, "FE", 26, 0},
      {43, "CDmso", 6, 0}, {17, "CH4", 6, 4}, {29, "CA2+", 20, 0},
      {22, "DUM", 0, 0}, {42, "SDmso", 16, 0}, {48, "CTFE", 6, 0},
      {5, "OW", 8, 0}, {34, "BR", 35, 0}, {49, "CHTFE", 6, 0}, {14, "CH1", 6, 1},
      {41, "HChl", 1, 0}, {53, "NUrea", 7, 0}, {23, "S", 16, 0},
      {13, "CH0", 6, 0}, {19, "CR1", 6, 1}, {11, "NE", 7, 0}, {16, "CH3", 6, 3},
      {8, "NL", 7, 0}, {12, "C", 6, 0}, {10, "NZ", 7, 0}, {6, "N", 7, 0},
      {24, "CU1+", 29, 0}, {50, "OTFE", 8, 0}, {25, "CU2+", 29, 0},
      {3, "OA", 8, 0}, {47, "FTFE", 9, 0}, {18, "CH2r", 6, 2}, {9, "NR", 7, 0}};
    for (auto& in : atom_dat)
      ff.NewAtomType(in.first, in.second, PT[in.third], in.fourth);
    
    // Add bond types
    std::vector<stdx::quad<int, double, double, double>> bnd_dat = {
      {34, 0.198, 50181.12, 640000.0}, {48, 0.290283, 502214.75, 2980000.0},
      {22, 0.148, 251019.84, 5730000.0}, {49, 0.279388, 373115.59, 2390000.0},
      {50, 0.291189, 371384.73, 2190000.0}, {10, 0.133, 417460.4, 11800000.0},
      {23, 0.148, 334693.12, 7640000.0}, {4, 0.112, 928256.0, 37000000.0},
      {14, 0.138, 418968.0, 11000000.0}, {31, 0.178, 376405.92, 5940000.0},
      {27, 0.153, 334748.7, 7150000.0}, {35, 0.2, 50240.0, 628000.0},
      {6, 0.125, 418750.0, 13400000.0}, {5, 0.123, 502282.8, 16600000.0},
      {7, 0.132, 418176.0, 12000000.0}, {12, 0.134, 420170.4, 11700000.0},
      {39, 0.11, 292820.0, 12100000.0}, {19, 0.143, 376670.58, 9210000.0},
      {28, 0.161, 250915.28, 4840000.0}, {21, 0.147, 376428.78, 8710000.0},
      {37, 0.221, 52748.28, 540000.0}, {46, 0.163299, 464531.53, 8710000.0},
      {45, 0.135, 375435.0, 10300000.0}, {16, 0.139, 417333.6, 10800000.0},
      {52, 0.287407, 502224.92, 3040000.0}, {41, 0.153, 376416.72, 8040000.0},
      {13, 0.136, 377318.4, 10200000.0}, {24, 0.148, 376748.8, 8600000.0},
      {25, 0.15, 376650.0, 8370000.0}, {17, 0.14, 334768.0, 8540000.0},
      {2, 0.1, 374000.0, 18700000.0}, {33, 0.187, 251077.42, 3590000.0},
      {36, 0.204, 418656.96, 5030000.0}, {15, 0.139, 334639.72, 8660000.0},
      {20, 0.1435, 251225.45, 6100000.0}, {43, 0.176, 501811.2, 8100000.0},
      {40, 0.1758, 501907.59, 8120000.0}, {32, 0.183, 376416.36, 5620000.0},
      {1, 0.1, 314000.0, 15700000.0}, {3, 0.109, 292272.6, 12300000.0},
      {30, 0.178, 172360.96, 2720000.0}, {11, 0.134, 377076.0, 10500000.0},
      {9, 0.133, 375006.8, 10600000.0}, {47, 0.233839, 293088.43, 2680000.0},
      {26, 0.152, 250909.44, 5430000.0}, {8, 0.133, 313802.86, 8870000.0},
      {42, 0.193799, 371824.72, 4950000.0}, {18, 0.143, 334545.64, 8180000.0},
      {29, 0.163, 250811.36, 4720000.0}, {38, 0.1, 464000.0, 23200000.0},
      {44, 0.1265, 419258.95, 13100000.0}, {51, 0.2077, 342525.96, 3970000.0}};
    for (auto& dat : bnd_dat)
      ff.LinkBondTypes(ff.NewBondType(BondType::Harmonic, dat.first, dat.third,
                                      dat.second),
                       ff.NewBondType(BondType::Quartic, dat.first, dat.fourth,
                                      dat.second));
    
    // Add angle types
    std::vector<stdx::quad<int, double, double, double>> ang_dat = {
      {1, 90.0, 0.11550101, 380.0}, {3, 96.0, 0.12177061, 405.0},
      {14, 109.6, 0.12142334, 450.0}, {9, 109.5, 0.08638614, 320.0},
      {40, 155.0, 0.12112698, 2215.0}, {15, 111.0, 0.14048747, 530.0},
      {29, 120.0, 0.17801113, 780.0}, {34, 125.0, 0.076490216, 375.0},
      {24, 120.0, 0.10147593, 445.0}, {37, 126.0, 0.12744672, 640.0},
      {45, 97.4, 0.14024534, 469.0}, {28, 120.0, 0.15288018, 670.0},
      {20, 116.0, 0.11421859, 465.0}, {12, 109.5, 0.12157397, 450.0},
      {27, 120.0, 0.12774922, 560.0}, {36, 126.0, 0.11448735, 575.0},
      {43, 107.57, 0.13376523, 484.0}, {46, 106.75, 0.14026005, 503.0},
      {51, 110.3, 0.14017954, 524.0}, {5, 103.0, 0.12122177, 420.0},
      {8, 109.5, 0.076912479, 285.0}, {47, 108.53, 0.12108416, 443.0},
      {19, 115.0, 0.1524165, 610.0}, {41, 180.0, 0.072640156, 91350.0},
      {16, 113.0, 0.14045138, 545.0}, {4, 100.0, 0.14008261, 475.0},
      {2, 90.0, 0.12768574, 420.0}, {50, 109.5, 0.12103261, 448.0},
      {18, 115.0, 0.11488482, 460.0}, {17, 115.0, 0.01229657, 50.0},
      {32, 123.0, 0.088743846, 415.0}, {44, 111.3, 0.16689058, 632.0},
      {39, 132.0, 0.12775497, 760.0}, {13, 109.5, 0.14052124, 520.0},
      {31, 122.0, 0.15317431, 700.0}, {53, 117.2, 0.15305438, 636.0},
      {42, 109.5, 0.11724316, 434.0}, {21, 116.0, 0.15236094, 620.0},
      {30, 121.0, 0.15312732, 685.0}, {35, 125.0, 0.1531408, 750.0},
      {6, 104.0, 0.14028506, 490.0}, {48, 109.5, 0.1670474, 618.0},
      {52, 111.4, 0.14025677, 532.0}, {54, 121.4, 0.1529482, 690.0},
      {7, 108.0, 0.12788754, 465.0}, {22, 117.0, 0.15336019, 635.0},
      {11, 109.5, 0.11480708, 425.0}, {26, 120.0, 0.12089532, 530.0},
      {10, 109.5, 0.10262668, 380.0}, {33, 124.0, 0.15266919, 730.0},
      {25, 120.0, 0.11518373, 505.0}, {38, 126.0, 0.15336544, 770.0},
      {49, 107.6, 0.14008648, 507.0}, {23, 120.0, 0.088910434, 390.0}};
    for (auto& dat : ang_dat)
      ff.LinkAngleTypes(ff.NewAngleType(AngleType::Harmonic, dat.first,
                                          dat.third, dat.second),
                        ff.NewAngleType(AngleType::CosineHarmonic, dat.first,
                                          dat.fourth, dat.second));
    
    // Add improper types
    std::vector<stdx::triple<int, double, double>> imp_dat = {
      {4, 0.051, 180.0}, {2, 0.102, 35.26439}, {5, 0.102, -35.26439},
      {3, 0.204, 0.0}, {1, 0.051, 0.0}};
    for (auto& dat : imp_dat)
      ff.NewDihedralType(DihedralType::Improper, dat.first, dat.second,
                         dat.third);
    
    // Add proper types
    std::vector<stdx::quad<int, double, double, int>> prp_dat = {
      {19, 3.14, 0.0, 2}, {42, 3.5, 180.0, 2}, {17, 0.418, 0.0, 2},
      {2, 3.41, 180.0, 1}, {7, 2.79, 0.0, 1}, {24, 1.3, 0.0, 3},
      {37, 9.5, 0.0, 3}, {45, 0.4, 0.0, 6}, {15, 41.8, 180.0, 2},
      {16, 0.0, 0.0, 2}, {38, 0.0, 0.0, 4}, {41, 3.77, 0.0, 6},
      {44, 0.7, 180.0, 6}, {22, 1.05, 0.0, 3}, {6, 9.45, 180.0, 1},
      {10, 5.86, 180.0, 2}, {28, 3.65, 0.0, 3}, {13, 24.0, 180.0, 2},
      {1, 2.67, 180.0, 1}, {35, 7.69, 0.0, 3}, {11, 7.11, 180.0, 2},
      {30, 3.9, 0.0, 3}, {3, 4.97, 180.0, 1}, {20, 5.09, 0.0, 2},
      {23, 1.26, 0.0, 3}, {32, 4.69, 0.0, 3}, {26, 2.93, 0.0, 3},
      {25, 2.53, 0.0, 3}, {9, 1.53, 180.0, 2}, {12, 16.7, 180.0, 2},
      {27, 3.19, 0.0, 3}, {14, 33.5, 180.0, 2}, {8, 5.35, 0.0, 1},
      {33, 5.44, 0.0, 3}, {34, 5.92, 0.0, 3}, {5, 9.35, 180.0, 1},
      {29, 3.77, 0.0, 3}, {39, 1.0, 180.0, 6}, {4, 5.86, 180.0, 1},
      {31, 4.18, 0.0, 3}, {43, 2.8, 0.0, 3}, {36, 8.62, 0.0, 3},
      {40, 1.0, 0.0, 6}, {21, 16.7, 0.0, 2}, {18, 2.09, 0.0, 2}};
    for (auto& dat : prp_dat)
      ff.NewDihedralType(DihedralType::Proper, dat.first, dat.second,
                          dat.third, dat.fourth);
    
    return ff;
  }
  
  // ===========================================================================
  // == Forcefield and Components Testing ======================================
  // ===========================================================================
  
  /*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFDihedral serilisation", T, ixffdhd_serial) {
   using In = typename T::t1;
   using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
   test::FFDihedralTestFixture fixture;
   
   std::ostringstream os;
   {
   Out oar(os);
   check_nothrow(oar(fixture.proper.imp, fixture.ff));
   }
   
   Forcefield loaded_ff;
   std::istringstream is(os.str());
   {
   In iar(is);
   check_nothrow(iar(fixture.dhd.imp, loaded_ff));
   }
   
   check_eq(fixture.proper.get_id(), fixture.dhd.get_id());
   check_eq(fixture.proper.get_dat(), fixture.dhd.get_dat());
   check_eq(fixture.proper.get_mask(), fixture.dhd.get_mask());
   check_eq(fixture.proper.get_type(), fixture.dhd.get_type());
   check_eq(loaded_ff, fixture.dhd.GetForcefield());
   }
   DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffdhd_serial, ixserial<IXFFDihedral>);
   */
  
  /*  test_case("IXFFDihedral construction") {
   Forcefield ff = test::CreateGenericTestForcefield().imp;
   using F = test::TestFFDihedral;
   using T = DihedralType;
   // Empty construction
   check_nothrow(F t(T::Empty, 4, {}, ff));
   check_nothrow(F t(T::Empty, ff));
   check_throws_as(F t(T::Empty, 4, {1.0}, ff), std::range_error);
   
   // Proper construction
   check_nothrow(F t(T::Proper, 4, {1.,2.,3}, ff));
   check_nothrow(F t(T::Proper, ff));
   check_throws_as(F t(T::Proper, 4, {1.,2.}, ff), std::range_error);
   check_throws_as(F t(T::Proper, 4, {1.,2.,3.,4.}, ff), std::range_error);
   
   // Improper construction
   check_nothrow(F t(T::Improper, 4, {1.,2.}, ff));
   check_nothrow(F t(T::Improper, ff));
   check_throws_as(F t(T::Improper, 4, {1.}, ff), std::range_error);
   check_throws_as(F t(T::Improper, 4, {1.,2.,3.}, ff), std::range_error);
   }
   
   test_case_fixture(test::FFDihedralTestFixture, "IXFFDihedral getting checks") {
   // Empty
   check_throws_as(empty.GetPhaseShift(), std::runtime_error);
   check_throws_as(empty.GetForceConstant(), std::runtime_error);
   check_throws_as(empty.GetMultiplicity(), std::runtime_error);
   check_throws_as(empty.GetIdealAngle(), std::runtime_error);
   check_eq(DihedralType::Empty, empty.GetType());
   check_eq(0, empty.GetID());
   check_eq(ff, empty.GetForcefield());
   
   // Proper
   check_eq(approximately(180.0), proper.GetPhaseShift());
   check_eq(approximately(2.67), proper.GetForceConstant());
   check_eq(1, proper.GetMultiplicity());
   check_throws_as(proper.GetIdealAngle(), std::runtime_error);
   check_eq(DihedralType::Proper, proper.GetType());
   check_eq(1, proper.GetID());
   check_eq(ff, proper.GetForcefield());
   
   // Improper
   check_throws_as(improper.GetPhaseShift(), std::runtime_error);
   check_eq(approximately(0.102), improper.GetForceConstant());
   check_throws_as(improper.GetMultiplicity(), std::runtime_error);
   check_eq(approximately(35.26439), improper.GetIdealAngle());
   check_eq(DihedralType::Improper, improper.GetType());
   check_eq(2, improper.GetID());
   check_eq(ff, improper.GetForcefield());
   }
   
   test_case_fixture(test::FFDihedralTestFixture, "IXFFDihedral internals checks") {
   using AllowEnum = indigox::test::TestFFDihedral::AllowEnum;
   // Empty
   check_eq(DihedralType::Empty, empty.get_type());
   check_eq(0, empty.get_id());
   check_eq(approximately(0.0), empty.get_dat()[0]);
   check_eq(approximately(0.0), empty.get_dat()[1]);
   check_eq(approximately(0.0), empty.get_dat()[2]);
   check_false(empty.get_mask()[AllowEnum::Allow_IdealAngle]);
   check_false(empty.get_mask()[AllowEnum::Allow_PhaseShift]);
   check_false(empty.get_mask()[AllowEnum::Allow_Multiplicity]);
   check_false(empty.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(ff, empty.get_ff().lock());
   
   // Proper
   check_eq(DihedralType::Proper, proper.get_type());
   check_eq(1, proper.get_id());
   check_eq(approximately(180.0), proper.get_dat()[0]);
   check_eq(approximately(2.67), proper.get_dat()[1]);
   check_eq(approximately(1.0), proper.get_dat()[2]);
   check_false(proper.get_mask()[AllowEnum::Allow_IdealAngle]);
   check(proper.get_mask()[AllowEnum::Allow_PhaseShift]);
   check(proper.get_mask()[AllowEnum::Allow_Multiplicity]);
   check(proper.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(ff, proper.get_ff().lock());
   
   // Improper
   check_eq(DihedralType::Improper, improper.get_type());
   check_eq(2, improper.get_id());
   check_eq(approximately(35.26439), improper.get_dat()[0]);
   check_eq(approximately(0.102), improper.get_dat()[1]);
   check_eq(approximately(0.0), improper.get_dat()[2]);
   check(improper.get_mask()[AllowEnum::Allow_IdealAngle]);
   check_false(improper.get_mask()[AllowEnum::Allow_PhaseShift]);
   check_false(improper.get_mask()[AllowEnum::Allow_Multiplicity]);
   check(improper.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(ff, improper.get_ff().lock());
   }
   */
  
  /*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAngle serilisation", T, ixffang_serial) {
   using In = typename T::t1;
   using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
   test::FFAngleTestFixture fixture;
   
   std::ostringstream os;
   {
   Out oar(os);
   check_nothrow(oar(fixture.harmonic.imp, fixture.ff));
   }
   
   Forcefield loaded_ff;
   std::istringstream is(os.str());
   {
   In iar(is);
   check_nothrow(iar(fixture.ang.imp, loaded_ff));
   }
   
   check_eq(fixture.harmonic.get_id(), fixture.ang.get_id());
   check_eq(fixture.harmonic.get_dat(), fixture.ang.get_dat());
   check_eq(fixture.harmonic.get_mask(), fixture.ang.get_mask());
   check_eq(fixture.harmonic.get_type(), fixture.ang.get_type());
   check_eq(loaded_ff, fixture.ang.GetForcefield());
   }
   DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffang_serial, ixserial<IXFFAngle>);
   */
  
  /*  test_case("IXFFAngle construction") {
   Forcefield ff = test::CreateGenericTestForcefield().imp;
   using F = test::TestFFAngle;
   using T = AngleType;
   // Empty construction
   check_nothrow(F t(T::Empty, 4, {}, ff));
   check_nothrow(F t(T::Empty, ff));
   check_throws_as(F t(T::Empty, 4, {1.}, ff), std::range_error);
   
   // Harmonic construction
   check_nothrow(F t(T::Harmonic, 4, {1.,2.}, ff));
   check_nothrow(F t(T::Harmonic, ff));
   check_throws_as(F t(T::Harmonic, 4, {1.}, ff), std::range_error);
   check_throws_as(F t(T::Empty, 4, {1.,2.,3.}, ff), std::range_error);
   
   // CosineHarmonic construction
   check_nothrow(F t(T::CosineHarmonic, 4, {1.,2.}, ff));
   check_nothrow(F t(T::CosineHarmonic, ff));
   check_throws_as(F t(T::CosineHarmonic, 4, {1.}, ff), std::range_error);
   check_throws_as(F t(T::CosineHarmonic, 4, {1.,2.,3.}, ff), std::range_error);
   }
   
   test_case_fixture(test::FFAngleTestFixture, "IXFFAngle getting checks") {
   // Empty
   check_throws_as(empty.GetIdealAngle(), std::runtime_error);
   check_throws_as(empty.GetForceConstant(), std::runtime_error);
   check_eq(AngleType::Empty, empty.GetType());
   check_eq(0, empty.GetID());
   check_eq(ff, empty.GetForcefield());
   
   // Harmonic
   check_eq(approximately(90.), harmonic.GetIdealAngle());
   check_eq(approximately(0.11550101), harmonic.GetForceConstant());
   check_eq(AngleType::Harmonic, harmonic.GetType());
   check_eq(1, harmonic.GetID());
   check_eq(ff, harmonic.GetForcefield());
   check_eq(cosineharmonic.imp, harmonic.GetLinkedType());
   
   // CosineHarmonic
   check_eq(approximately(90.), cosineharmonic.GetIdealAngle());
   check_eq(approximately(380.), cosineharmonic.GetForceConstant());
   check_eq(AngleType::CosineHarmonic, cosineharmonic.GetType());
   check_eq(1, cosineharmonic.GetID());
   check_eq(ff, cosineharmonic.GetForcefield());
   check_eq(harmonic.imp, cosineharmonic.GetLinkedType());
   }
   
   test_case_fixture(test::FFAngleTestFixture, "IXFFAngle internals checks") {
   using AllowEnum = indigox::test::TestFFAngle::AllowEnum;
   // Empty
   check_eq(AngleType::Empty, empty.get_type());
   check_eq(0, empty.get_id());
   check_eq(approximately(0.0), empty.get_dat()[0]);
   check_eq(approximately(0.0), empty.get_dat()[1]);
   check_false(empty.get_mask()[AllowEnum::Allow_IdealAngle]);
   check_false(empty.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(FFAngle(), empty.get_link());
   check_eq(ff, empty.get_ff().lock());
   
   // Harmonic
   check_eq(AngleType::Harmonic, harmonic.get_type());
   check_eq(1, harmonic.get_id());
   check_eq(approximately(0.11550101), harmonic.get_dat()[0]);
   check_eq(approximately(90.), harmonic.get_dat()[1]);
   check(harmonic.get_mask()[AllowEnum::Allow_IdealAngle]);
   check(harmonic.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(cosineharmonic.imp, harmonic.get_link());
   check_eq(ff, harmonic.get_ff().lock());
   
   // CosineHarmonic
   check_eq(AngleType::CosineHarmonic, cosineharmonic.get_type());
   check_eq(1, cosineharmonic.get_id());
   check_eq(approximately(380.), cosineharmonic.get_dat()[0]);
   check_eq(approximately(90.), cosineharmonic.get_dat()[1]);
   check(cosineharmonic.get_mask()[AllowEnum::Allow_IdealAngle]);
   check(cosineharmonic.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(harmonic.imp, cosineharmonic.get_link());
   check_eq(ff, cosineharmonic.get_ff().lock());
   }
   */
  
  /*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFBond serilisation", T, ixffbnd_serial) {
   using In = typename T::t1;
   using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
   test::FFBondTestFixture fixture;
   
   std::ostringstream os;
   {
   Out oar(os);
   check_nothrow(oar(fixture.harmonic.imp, fixture.ff));
   }
   
   Forcefield loaded_ff;
   std::istringstream is(os.str());
   {
   In iar(is);
   check_nothrow(iar(fixture.bnd.imp, loaded_ff));
   }
   
   check_eq(fixture.harmonic.get_id(), fixture.bnd.get_id());
   check_eq(fixture.harmonic.get_dat(), fixture.bnd.get_dat());
   check_eq(fixture.harmonic.get_mask(), fixture.bnd.get_mask());
   check_eq(fixture.harmonic.get_type(), fixture.bnd.get_type());
   check_eq(loaded_ff, fixture.bnd.GetForcefield());
   }
   DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffbnd_serial, ixserial<IXFFBond>);
   */
  
  /*  test_case("IXFFBond construction") {
   Forcefield ff = test::CreateGenericTestForcefield().imp;
   using F = test::TestFFBond;
   using T = BondType;
   // Empty construction
   check_nothrow(F t(T::Empty, 4, {}, ff));
   check_nothrow(F t(T::Empty, ff));
   check_throws_as(F t(T::Empty, 4, {1.}, ff), std::range_error);
   
   // Harmonic construction
   check_nothrow(F t(T::Harmonic, 4, {1.,2.}, ff));
   check_nothrow(F t(T::Harmonic, ff));
   check_throws_as(F t(T::Harmonic, 4, {1.}, ff), std::range_error);
   check_throws_as(F t(T::Empty, 4, {1.,2.,3.}, ff), std::range_error);
   
   // Quartic construction
   check_nothrow(F t(T::Quartic, 4, {1.,2.}, ff));
   check_nothrow(F t(T::Quartic, ff));
   check_throws_as(F t(T::Quartic, 4, {1.}, ff), std::range_error);
   check_throws_as(F t(T::Quartic, 4, {1.,2.,3.}, ff), std::range_error);
   }
   
   test_case_fixture(test::FFBondTestFixture, "IXFFBond getting checks") {
   // Empty
   check_throws_as(empty.GetIdealLength(), std::runtime_error);
   check_throws_as(empty.GetForceConstant(), std::runtime_error);
   check_eq(BondType::Empty, empty.GetType());
   check_eq(0, empty.GetID());
   check_eq(ff, empty.GetForcefield());
   
   // Harmonic
   check_eq(approximately(0.109), harmonic.GetIdealLength());
   check_eq(approximately(292272.6), harmonic.GetForceConstant());
   check_eq(BondType::Harmonic, harmonic.GetType());
   check_eq(3, harmonic.GetID());
   check_eq(ff, harmonic.GetForcefield());
   check_eq(quartic.imp, harmonic.GetLinkedType());
   
   // Quartic
   check_eq(approximately(0.109), quartic.GetIdealLength());
   check_eq(approximately(12300000.), quartic.GetForceConstant());
   check_eq(BondType::Quartic, quartic.GetType());
   check_eq(3, quartic.GetID());
   check_eq(ff, quartic.GetForcefield());
   check_eq(harmonic.imp, quartic.GetLinkedType());
   }
   
   test_case_fixture(test::FFBondTestFixture, "IXFFBond internals checks") {
   using AllowEnum = indigox::test::TestFFBond::AllowEnum;
   // Empty
   check_eq(BondType::Empty, empty.get_type());
   check_eq(0, empty.get_id());
   check_eq(approximately(0.0), empty.get_dat()[0]);
   check_eq(approximately(0.0), empty.get_dat()[1]);
   check_false(empty.get_mask()[AllowEnum::Allow_IdealLength]);
   check_false(empty.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(FFBond(), empty.get_link());
   check_eq(ff, empty.get_ff().lock());
   
   // Harmonic
   check_eq(BondType::Harmonic, harmonic.get_type());
   check_eq(3, harmonic.get_id());
   check_eq(approximately(292272.6), harmonic.get_dat()[0]);
   check_eq(approximately(0.109), harmonic.get_dat()[1]);
   check(harmonic.get_mask()[AllowEnum::Allow_IdealLength]);
   check(harmonic.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(quartic.imp, harmonic.get_link());
   check_eq(ff, harmonic.get_ff().lock());
   
   // Quartic
   check_eq(BondType::Quartic, quartic.get_type());
   check_eq(3, quartic.get_id());
   check_eq(approximately(12300000.), quartic.get_dat()[0]);
   check_eq(approximately(0.109), quartic.get_dat()[1]);
   check(quartic.get_mask()[AllowEnum::Allow_IdealLength]);
   check(quartic.get_mask()[AllowEnum::Allow_ForceConstant]);
   check_eq(harmonic.imp, quartic.get_link());
   check_eq(ff, quartic.get_ff().lock());
   }
   */
  
  /*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("FFAtom serialisation", T, ixffatm_serial) {
   using In = typename T::t1;
   using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
   test::FFAtomTestFixture fixture;
   
   std::ostringstream os;
   {
   Out oar(os);
   check_nothrow(oar(fixture.atm.imp, fixture.ff));
   }
   
   test::TestFFAtom loaded(0, "", fixture.ff);
   Forcefield loaded_ff;
   std::istringstream is(os.str());
   {
   In iar(is);
   check_nothrow(iar(loaded.imp, loaded_ff));
   }
   
   check_eq(fixture.atm.get_id(), loaded.get_id());
   check_eq(fixture.atm.get_name(), loaded.get_name());
   check_eq(loaded_ff, loaded.get_ff().lock());
   }
   DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixffatm_serial, ixserial<IXFFAtom>);
   */
  
  /*  test_case("IXFFAtom construction") {
   Forcefield ff = test::CreateGenericTestForcefield().imp;
   using F = test::TestFFAtom;
   check_nothrow(F t(23, "TEST", ff));
   }
   
   test_case_fixture(test::FFAtomTestFixture, "IXFFAtom getting checks") {
   check_eq(7, atm.GetID());
   check_eq("ATM", atm.GetName());
   check_eq(ff, atm.GetForcefield());
   }
   
   test_case_fixture(test::FFAtomTestFixture, "IXFFAtom internals checks") {
   check_eq(7, atm.get_id());
   check_eq("ATM", atm.get_name());
   check_eq(ff, atm.get_ff().lock());
   }
   */
}
