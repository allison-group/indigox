#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/quad.hpp>
#include <indigox/utils/serialise.hpp>
#include <indigox/utils/triple.hpp>

#include <array>
#include <initializer_list>
#include <numeric>

namespace indigox {

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid ffatom instance")
#else
#define _sanity_check_(x)
#endif

  // ===========================================================================
  // == FFAtom Data Implementation =============================================
  // ===========================================================================

  struct FFAtom::Impl {
    int32_t id;
    std::string name;
    int32_t element;
    Forcefield ff;

    Impl() = default;
    Impl(int32_t i, std::string n, const Element &e, const Forcefield &f)
        : id(i), name(n), element(e.GetAtomicNumber()), ff(f) {}

    template <class Archive> void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("id", id), INDIGOX_SERIAL_NVP("name", name),
              INDIGOX_SERIAL_NVP("element", element),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };

  // ===========================================================================
  // == FFAtom Construction and Assignment =====================================
  // ===========================================================================

  FFAtom::FFAtom(int32_t id, std::string name, const Element &e,
                 const Forcefield &ff)
      : m_data(std::make_shared<Impl>(id, name, e, ff)) {}

  // ===========================================================================
  // == FFAtom Serialisation ===================================================
  // ===========================================================================

  template <class Archive>
  void FFAtom::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(FFAtom);

  // ===========================================================================
  // == FFAtom Data Retrevial ==================================================
  // ===========================================================================

  int32_t FFAtom::GetID() const {
    _sanity_check_(*this);
    return m_data->id;
  }
  std::string FFAtom::GetName() const {
    _sanity_check_(*this);
    return m_data->name;
  }
  Forcefield FFAtom::GetForcefield() const {
    _sanity_check_(*this);
    return m_data->ff;
  }
  Element FFAtom::GetElement() const {
    _sanity_check_(*this);
    return GetPeriodicTable()[m_data->element];
  }

  // ===========================================================================
  // == FFAtom Operators =======================================================
  // ===========================================================================

  bool FFAtom::operator==(const FFAtom &atm) const {
    _sanity_check_(*this);
    _sanity_check_(atm);
    return GetID() == atm.GetID() && GetName() == atm.GetName();
  }

  bool FFAtom::operator<(const FFAtom &atm) const {
    _sanity_check_(*this);
    _sanity_check_(atm);
    return GetID() < atm.GetID();
  }

  bool FFAtom::operator>(const FFAtom &atm) const {
    _sanity_check_(*this);
    _sanity_check_(atm);
    return GetID() > atm.GetID();
  }

  std::ostream &operator<<(std::ostream &os, const FFAtom &atm) {
    if (atm) { os << "FFAtom(" << atm.GetID() << ")"; }
    return os;
  }

#undef _sanity_check_

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid ffbond instance")
#else
#define _sanity_check_(x)
#endif

  // ===========================================================================
  // == FFBond Data implementation =============================================
  // ===========================================================================

  struct FFBond::Impl {
    BondType type;
    int32_t id;
    FFBond::DataStore raw_data;
    FFBond::AllowedMask allowed_parameters;
    FFBond linked_bond;
    Forcefield ff;

    Impl() = default;
    Impl(BondType t, int32_t i, FFParam data, const Forcefield &f)
        : type(t), id(i), linked_bond(), ff(f) {
      if (t == BondType::Harmonic) allowed_parameters.from_uint32(3);
      if (t == BondType::Quartic) allowed_parameters.from_uint32(3);
      raw_data.fill(0);
      if (data.size() != allowed_parameters.count())
        throw std::range_error("Incorrect item count");
      std::copy(data.begin(), data.end(), raw_data.begin());
    }

    template <class Archive> void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("type", type), INDIGOX_SERIAL_NVP("id", id),
              INDIGOX_SERIAL_NVP("raw_data", raw_data),
              INDIGOX_SERIAL_NVP("allowed", allowed_parameters),
              INDIGOX_SERIAL_NVP("linked_bond_type", linked_bond),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };

  // ===========================================================================
  // == FFBond Construction and Assignment =====================================
  // ===========================================================================

  FFBond::FFBond(BondType type, int32_t id, FFParam params,
                 const Forcefield &ff)
      : m_data(std::make_shared<Impl>(type, id, params, ff)) {}

  // ===========================================================================
  // == FFBond Serialisation ===================================================
  // ===========================================================================

  template <class Archive>
  void FFBond::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(FFBond)

  // ===========================================================================
  // == FFBond Data Retrevial ==================================================
  // ===========================================================================

  double FFBond::GetForceConstant() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_ForceConstant]
            ? m_data->raw_data[Store_ForceConstant]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFBond::GetIdealLength() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_IdealLength]
            ? m_data->raw_data[Store_IdealLength]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  BondType FFBond::GetType() const {
    _sanity_check_(*this);
    return m_data->type;
  }
  int32_t FFBond::GetID() const {
    _sanity_check_(*this);
    return m_data->id;
  }
  FFBond FFBond::GetLinkedType() const {
    _sanity_check_(*this);
    return m_data->linked_bond;
  }
  Forcefield FFBond::GetForcefield() const {
    _sanity_check_(*this);
    return m_data->ff;
  }

  // ===========================================================================
  // == FFBond Operators =======================================================
  // ===========================================================================

  bool FFBond::operator==(const FFBond &bnd) const {
    _sanity_check_(*this);
    _sanity_check_(bnd);
    return (GetType() == bnd.GetType() && GetID() == bnd.GetID() &&
            m_data->raw_data == bnd.m_data->raw_data);
  }

  bool FFBond::operator<(const FFBond &bnd) const {
    _sanity_check_(*this);
    _sanity_check_(bnd);
    if (GetType() < bnd.GetType()) return true;
    if (GetID() < bnd.GetID()) return true;
    return false;
  }
  bool FFBond::operator>(const FFBond &bnd) const {
    _sanity_check_(*this);
    _sanity_check_(bnd);
    if (GetType() > bnd.GetType()) return true;
    if (GetID() > bnd.GetID()) return true;
    return false;
  }

  std::ostream &operator<<(std::ostream &os, BondType type) {
    switch (type) {
    case FFBond::Type::Empty: return os << "Empty";
    case FFBond::Type::Harmonic: return os << "Harmonic";
    case FFBond::Type::Quartic: return os << "Quartic";
    default: return os;
    }
  }

  std::ostream &operator<<(std::ostream &os, const FFBond &bnd) {
    if (bnd) os << "FFBond(" << bnd.GetType() << ", " << bnd.GetID() << ")";
    return os;
  }

#undef _sanity_check_

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid ffangle instance")
#else
#define _sanity_check_(x)
#endif

  // ===========================================================================
  // == FFAngle Data Implementation ============================================
  // ===========================================================================

  struct FFAngle::Impl {
    AngleType type;
    int32_t id;
    FFAngle::DataStore raw_data;
    FFAngle::AllowedMask allowed_parameters;
    FFAngle linked_angle;
    Forcefield ff;

    Impl() = default;
    Impl(AngleType t, int32_t i, FFParam data, const Forcefield &f)
        : type(t), id(i), linked_angle(), ff(f) {
      if (t == AngleType::Harmonic) allowed_parameters.from_uint32(3);
      if (t == AngleType::CosineHarmonic) allowed_parameters.from_uint32(3);
      raw_data.fill(0);
      if (data.size() != allowed_parameters.count())
        throw std::range_error("Incorrect item count");
      std::copy(data.begin(), data.end(), raw_data.begin());
    }

    template <class Archive> void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("type", type), INDIGOX_SERIAL_NVP("id", id),
              INDIGOX_SERIAL_NVP("raw_data", raw_data),
              INDIGOX_SERIAL_NVP("allwed_mask", allowed_parameters),
              INDIGOX_SERIAL_NVP("link", linked_angle),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };

  // ===========================================================================
  // == FFAngle Construction and Assignment ====================================
  // ===========================================================================

  FFAngle::FFAngle(AngleType type, int32_t id, FFParam par,
                   const Forcefield &ff)
      : m_data(std::make_shared<Impl>(type, id, par, ff)) {}

  // ===========================================================================
  // == FFAngle Serialisation ==================================================
  // ===========================================================================

  template <class Archive>
  void FFAngle::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(FFAngle);

  // ===========================================================================
  // == FFAngle Data Retrevial =================================================
  // ===========================================================================

  double FFAngle::GetForceConstant() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_ForceConstant]
            ? m_data->raw_data[Store_ForceConstant]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFAngle::GetIdealAngle() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_IdealAngle]
            ? m_data->raw_data[Store_IdealAngle]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  AngleType FFAngle::GetType() const {
    _sanity_check_(*this);
    return m_data->type;
  }
  int32_t FFAngle::GetID() const {
    _sanity_check_(*this);
    return m_data->id;
  }
  FFAngle FFAngle::GetLinkedType() const {
    _sanity_check_(*this);
    return m_data->linked_angle;
  }
  Forcefield FFAngle::GetForcefield() const {
    _sanity_check_(*this);
    return m_data->ff;
  }

  // ===========================================================================
  // == FFAngle Operators ======================================================
  // ===========================================================================

  bool FFAngle::operator==(const FFAngle &ang) const {
    _sanity_check_(*this);
    _sanity_check_(ang);
    return (GetType() == ang.GetType() && GetID() == ang.GetID() &&
            m_data->raw_data == ang.m_data->raw_data);
  }
  bool FFAngle::operator<(const FFAngle &ang) const {
    _sanity_check_(*this);
    _sanity_check_(ang);
    if (GetType() < ang.GetType()) return true;
    if (GetID() < ang.GetID()) return true;
    return false;
  }
  bool FFAngle::operator>(const FFAngle &ang) const {
    _sanity_check_(*this);
    _sanity_check_(ang);
    if (GetType() > ang.GetType()) return true;
    if (GetID() > ang.GetID()) return true;
    return false;
  }

  std::ostream &operator<<(std::ostream &os, AngleType type) {
    switch (type) {
    case FFAngle::Type::Empty: return os << "Empty";
    case FFAngle::Type::Harmonic: return os << "Harmonic";
    case FFAngle::Type::CosineHarmonic: return os << "CosineHarmonic";
    default: return os;
    }
  }

  std::ostream &operator<<(std::ostream &os, const FFAngle &ang) {
    if (ang) os << "FFAngle(" << ang.GetType() << ", " << ang.GetID() << ")";
    return os;
  }

#undef _sanity_check_

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid ffdihedral instance")
#else
#define _sanity_check_(x)
#endif

  // ===========================================================================
  // == FFDihedral Data Implementation =========================================
  // ===========================================================================

  struct FFDihedral::Impl {
    DihedralType type;
    int32_t id;
    FFDihedral::DataStore raw_data;
    FFDihedral::AllowedMask allowed_parameters;
    Forcefield ff;

    Impl() = default;
    Impl(DihedralType t, int32_t i, FFParam data, const Forcefield &f)
        : type(t), id(i), ff(f) {
      if (t == DihedralType::Proper) allowed_parameters.from_uint32(7);
      if (t == DihedralType::Improper) allowed_parameters.from_uint32(10);
      raw_data.fill(0);
      if (data.size() != allowed_parameters.count())
        throw std::range_error("Incorrect item count");
      std::copy(data.begin(), data.end(), raw_data.begin());
    }

    template <class Archive> void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("type", type), INDIGOX_SERIAL_NVP("id", id),
              INDIGOX_SERIAL_NVP("raw_data", raw_data),
              INDIGOX_SERIAL_NVP("allwed_mask", allowed_parameters),
              INDIGOX_SERIAL_NVP("ff", ff));
    }
  };

  // ===========================================================================
  // == FFDihedral Construction and Assignment
  // ====================================
  // ===========================================================================

  FFDihedral::FFDihedral(DihedralType type, int32_t id, FFParam par,
                         const Forcefield &ff)
      : m_data(std::make_shared<Impl>(type, id, par, ff)) {}

  // ===========================================================================
  // == FFDihedral Serialisation ===============================================
  // ===========================================================================

  template <class Archive>
  void FFDihedral::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(FFDihedral);

  // ===========================================================================
  // == FFDihedral Data Retrevial ==============================================
  // ===========================================================================

  double FFDihedral::GetPhaseShift() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_PhaseShift]
            ? m_data->raw_data[Store_PhaseShift]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFDihedral::GetForceConstant() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_ForceConstant]
            ? m_data->raw_data[Store_ForceConstant]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  int32_t FFDihedral::GetMultiplicity() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_Multiplicity]
            ? static_cast<int32_t>(m_data->raw_data[Store_Multiplicity])
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  double FFDihedral::GetIdealAngle() const {
    _sanity_check_(*this);
    return (
        m_data->allowed_parameters[Allow_IdealAngle]
            ? m_data->raw_data[Store_IdealAngle]
            : throw std::runtime_error("Disallowed parameter type requested"));
  }
  DihedralType FFDihedral::GetType() const {
    _sanity_check_(*this);
    return m_data->type;
  }
  int32_t FFDihedral::GetID() const {
    _sanity_check_(*this);
    return m_data->id;
  }
  Forcefield FFDihedral::GetForcefield() const {
    _sanity_check_(*this);
    return m_data->ff;
  }

  // ===========================================================================
  // == FFDihedral Operators ===================================================
  // ===========================================================================

  bool FFDihedral::operator==(const FFDihedral &dhd) const {
    _sanity_check_(*this);
    _sanity_check_(dhd);
    return (GetType() == dhd.GetType() && GetID() == dhd.GetID() &&
            m_data->raw_data == dhd.m_data->raw_data);
  }
  bool FFDihedral::operator<(const FFDihedral &dhd) const {
    _sanity_check_(*this);
    _sanity_check_(dhd);
    if (GetType() < dhd.GetType()) return true;
    if (GetID() < dhd.GetID()) return true;
    return false;
  }
  bool FFDihedral::operator>(const FFDihedral &dhd) const {
    _sanity_check_(*this);
    _sanity_check_(dhd);
    if (GetType() > dhd.GetType()) return true;
    if (GetID() > dhd.GetID()) return true;
    return false;
  }

  std::ostream &operator<<(std::ostream &os, DihedralType dhd) {
    switch (dhd) {
    case FFDihedral::Type::Empty: return (os << "Empty");
    case FFDihedral::Type::Proper: return (os << "Proper");
    case FFDihedral::Type::Improper: return (os << "Improper");
    default: return os;
    }
  }

  std::ostream &operator<<(std::ostream &os, const FFDihedral &dhd) {
    if (dhd) os << "FFDihedral(" << dhd.GetType() << ", " << dhd.GetID() << ")";
    return os;
  }

#undef _sanity_check_
#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid forcefield instance")
#else
#define _sanity_check_(x)
#endif

  // ===========================================================================
  // == Forcefield Data Implementation =========================================
  // ===========================================================================

  struct Forcefield::Impl {
    FFFamily family;
    std::string name;
    Forcefield::AtomTypes atoms;
    Forcefield::BondTypes bonds;
    Forcefield::AngleTypes angles;
    Forcefield::DihedralTypes dihedrals;

    Impl() = default;
    Impl(FFFamily fam, std::string nme) : family(fam), name(nme) {
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

    template <class Archive> void serialise(Archive &archive, const uint32_t) {
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
                           [&i](auto &b) { return b.GetID() == i; }) !=
              bonds.at(t).end());
    }
    bool Contains(AngleType t) { return angles.find(t) != angles.end(); }
    bool Contains(AngleType t, int32_t i) {
      if (!Contains(t)) return false;
      return (std::find_if(angles.at(t).begin(), angles.at(t).end(),
                           [&i](auto &b) { return b.GetID() == i; }) !=
              angles.at(t).end());
    }
    bool Contains(DihedralType t) {
      return dihedrals.find(t) != dihedrals.end();
    }
    bool Contains(DihedralType t, int32_t i) {
      if (!Contains(t)) return false;
      return (std::find_if(dihedrals.at(t).begin(), dihedrals.at(t).end(),
                           [&i](auto &b) { return b.GetID() == i; }) !=
              dihedrals.at(t).end());
    }
  };

  // ===========================================================================
  // == Forcefield Construction and Assignment =================================
  // ===========================================================================

  Forcefield::Forcefield(FFFamily family, std::string name)
      : m_data(std::make_shared<Impl>(family, name)) {}

  // ===========================================================================
  // == Forcefield Serialisation ===============================================
  // ===========================================================================

  template <class Archive>
  void Forcefield::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(Forcefield);

  // ===========================================================================
  // == Forcefield Data Modification ===========================================
  // ===========================================================================

  FFBond Forcefield::NewBondType(BondType type, int32_t id, FFParam param) {
    _sanity_check_(*this);
    if (!m_data->Contains(type))
      throw std::out_of_range("Unsupported bond type");
    if (m_data->Contains(type, id))
      throw std::out_of_range("Bond ID already exists");
    m_data->bonds[type].emplace_back(FFBond(type, id, param, *this));
    return m_data->bonds[type].back();
  }

  FFAngle Forcefield::NewAngleType(AngleType type, int32_t id, FFParam param) {
    _sanity_check_(*this);
    if (!m_data->Contains(type))
      throw std::out_of_range("Unsupported angle type");
    if (m_data->Contains(type, id))
      throw std::out_of_range("Angle ID already exists");
    m_data->angles[type].emplace_back(FFAngle(type, id, param, *this));
    return m_data->angles[type].back();
  }

  FFDihedral Forcefield::NewDihedralType(DihedralType type, int32_t id,
                                         FFParam param) {
    _sanity_check_(*this);
    if (!m_data->Contains(type))
      throw std::out_of_range("Unsupported dihedral type");
    if (m_data->Contains(type, id))
      throw std::out_of_range("Dihedral ID already exists");
    m_data->dihedrals[type].emplace_back(FFDihedral(type, id, param, *this));
    return m_data->dihedrals[type].back();
  }

  FFAtom Forcefield::NewAtomType(int32_t id, std::string name,
                                 const Element &element) {
    _sanity_check_(*this);
    FFAtom atm(id, name, element, *this);
    if (std::find(m_data->atoms.begin(), m_data->atoms.end(), atm) !=
        m_data->atoms.end())
      throw std::out_of_range("Atom type already exists");
    m_data->atoms.emplace_back(atm);
    return m_data->atoms.back();
  }

  void Forcefield::ReserveAtomTypes(size_t sz) { m_data->atoms.reserve(sz); }

  FFBond Forcefield::NewBondType(BondType t, int32_t i, double a, double b) {
    return NewBondType(t, i, {a, b});
  }

  void Forcefield::LinkBondTypes(const FFBond &a, const FFBond &b) {
    _sanity_check_(*this);
    _sanity_check_(a);
    _sanity_check_(b);
    a.m_data->linked_bond = b;
    b.m_data->linked_bond = a;
  }

  void Forcefield::ReserveBondTypes(BondType type, size_t size) {
    if (m_data->Contains(type)) m_data->bonds.at(type).reserve(size);
  }

  FFAngle Forcefield::NewAngleType(AngleType t, int32_t i, double a, double b) {
    return NewAngleType(t, i, {a, b});
  }

  void Forcefield::LinkAngleTypes(const FFAngle &a, const FFAngle &b) {
    _sanity_check_(*this);
    _sanity_check_(a);
    _sanity_check_(b);
    a.m_data->linked_angle = b;
    b.m_data->linked_angle = a;
  }

  void Forcefield::ReserveAngleTypes(AngleType type, size_t size) {
    if (m_data->Contains(type)) m_data->angles.at(type).reserve(size);
  }

  FFDihedral Forcefield::NewDihedralType(DihedralType type, int32_t id,
                                         double a, double b, double c) {
    return NewDihedralType(type, id, {b, a, c});
  }

  FFDihedral Forcefield::NewDihedralType(DihedralType type, int32_t id,
                                         double a, double b) {
    return NewDihedralType(type, id, {b, a});
  }

  void Forcefield::ReserveDihedralTypes(DihedralType type, size_t size) {
    _sanity_check_(*this);
    if (m_data->Contains(type)) m_data->dihedrals.at(type).reserve(size);
  }

  // ================================================================
  // == Forcefield Data Retrevial ===================================
  // ================================================================

  template <typename InputIt> size_t sum_of_sizes(InputIt begin, InputIt end) {
    return std::accumulate(
        begin, end, 0, [](size_t a, auto &b) { return a + b.second.size(); });
  }

  FFAtom Forcefield::GetAtomType(std::string name) const {
    _sanity_check_(*this);
    auto bgn = m_data->atoms.begin();
    auto end = m_data->atoms.end();
    auto fnd = [&name](FFAtom atm) { return atm.GetName() == name; };
    auto pos = std::find_if(bgn, end, fnd);
    return (pos == end) ? FFAtom() : *pos;
  }

  FFAtom Forcefield::GetAtomType(int32_t id) const {
    _sanity_check_(*this);
    auto bgn = m_data->atoms.begin();
    auto end = m_data->atoms.end();
    auto fnd = [&id](auto &Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    return (pos == end) ? FFAtom() : *pos;
  }

  size_t Forcefield::NumAtomTypes() const { return m_data->atoms.size(); }

  FFBond Forcefield::GetBondType(BondType type, int32_t id) const {
    _sanity_check_(*this);
    if (!m_data->Contains(type)) { return FFBond(); }
    auto bgn = m_data->bonds.at(type).begin();
    auto end = m_data->bonds.at(type).end();
    auto fnd = [&id](auto &Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    return (pos == end) ? FFBond() : *pos;
  }
  FFBond Forcefield::GetBondType(int32_t id) const {
    _sanity_check_(*this);
    FFBond typ;
    for (auto &types : m_data->bonds) {
      typ = GetBondType(types.first, id);
      if (typ) { break; }
    }
    return typ;
  }

  size_t Forcefield::NumBondTypes() const {
    _sanity_check_(*this);
    return sum_of_sizes(m_data->bonds.begin(), m_data->bonds.end());
  }

  size_t Forcefield::NumBondTypes(BondType type) const {
    _sanity_check_(*this);
    return m_data->Contains(type) ? m_data->bonds.at(type).size() : 0;
  }
  FFAngle Forcefield::GetAngleType(AngleType type, int32_t id) const {
    _sanity_check_(*this);
    if (!m_data->Contains(type)) { return FFAngle(); }
    auto bgn = m_data->angles.at(type).begin();
    auto end = m_data->angles.at(type).end();
    auto fnd = [&id](auto &Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    return (pos == end) ? FFAngle() : *pos;
  }

  FFAngle Forcefield::GetAngleType(int id) const {
    _sanity_check_(*this);
    FFAngle ang;
    for (auto &types : m_data->angles) {
      ang = GetAngleType(types.first, id);
      if (ang) { break; }
    }
    return ang;
  }

  size_t Forcefield::NumAngleTypes() const {
    _sanity_check_(*this);
    return sum_of_sizes(m_data->angles.begin(), m_data->angles.end());
  }

  size_t Forcefield::NumAngleTypes(AngleType type) const {
    _sanity_check_(*this);
    return m_data->Contains(type) ? m_data->angles.at(type).size() : 0;
  }

  FFDihedral Forcefield::GetDihedralType(DihedralType type, int32_t id) const {
    _sanity_check_(*this);
    if (!m_data->Contains(type)) { return FFDihedral(); }
    auto bgn = m_data->dihedrals.at(type).begin();
    auto end = m_data->dihedrals.at(type).end();
    auto fnd = [&id](auto &Z) { return Z.GetID() == id; };
    auto pos = std::find_if(bgn, end, fnd);
    return (pos == end) ? FFDihedral() : *pos;
  }

  FFDihedral Forcefield::GetDihedralType(int32_t id) const {
    _sanity_check_(*this);
    FFDihedral dhd;
    for (auto &types : m_data->dihedrals) {
      dhd = GetDihedralType(types.first, id);
      if (dhd) { break; }
    }
    return dhd;
  }

  size_t Forcefield::NumDihedralTypes() const {
    _sanity_check_(*this);
    return sum_of_sizes(m_data->dihedrals.begin(), m_data->dihedrals.end());
  }

  size_t Forcefield::NumDihedralTypes(DihedralType type) const {
    _sanity_check_(*this);
    return m_data->Contains(type) ? m_data->dihedrals.at(type).size() : 0;
  }
  FFFamily Forcefield::GetFamily() const {
    _sanity_check_(*this);
    return m_data->family;
  }
  std::string Forcefield::GetName() const {
    _sanity_check_(*this);
    return m_data->name;
  }

  // ===========================================================================
  // == Forcefield Operators ===================================================
  // ===========================================================================

  bool Forcefield::operator==(const Forcefield &ff) const {
    _sanity_check_(*this);
    _sanity_check_(ff);
    if (GetFamily() != ff.GetFamily()) return false;
    if (GetName() != ff.GetName()) return false;
    /// \todo Add more indepth comparision
    return true;
  }

  std::ostream &operator<<(std::ostream &os, const Forcefield &ff) {
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
        {2, "OM", 8, 0},      {1, "O", 8, 0},       {4, "OE", 8, 0},
        {39, "CChl", 6, 0},   {32, "F", 9, 0},      {52, "OUrea", 8, 0},
        {15, "CH2", 6, 2},    {20, "HC", 1, 0},     {37, "NA+", 11, 0},
        {33, "CL", 17, 0},    {44, "ODmso", 8, 0},  {40, "CLChl", 17, 0},
        {45, "CCl4", 6, 0},   {31, "AR", 18, 0},    {38, "CL-", 17, 0},
        {35, "CMet", 6, 0},   {21, "H", 1, 0},      {27, "ZN2+", 30, 0},
        {28, "MG2+", 12, 0},  {46, "CLCl4", 17, 0}, {36, "OMet", 8, 0},
        {54, "CH3p", 6, 3},   {51, "CUrea", 6, 0},  {30, "P,SI", 15, 0},
        {7, "NT", 7, 0},      {26, "FE", 26, 0},    {43, "CDmso", 6, 0},
        {17, "CH4", 6, 4},    {29, "CA2+", 20, 0},  {22, "DUM", 0, 0},
        {42, "SDmso", 16, 0}, {48, "CTFE", 6, 0},   {5, "OW", 8, 0},
        {34, "BR", 35, 0},    {49, "CHTFE", 6, 0},  {14, "CH1", 6, 1},
        {41, "HChl", 1, 0},   {53, "NUrea", 7, 0},  {23, "S", 16, 0},
        {13, "CH0", 6, 0},    {19, "CR1", 6, 1},    {11, "NE", 7, 0},
        {16, "CH3", 6, 3},    {8, "NL", 7, 0},      {12, "C", 6, 0},
        {10, "NZ", 7, 0},     {6, "N", 7, 0},       {24, "CU1+", 29, 0},
        {50, "OTFE", 8, 0},   {25, "CU2+", 29, 0},  {3, "OA", 8, 0},
        {47, "FTFE", 9, 0},   {18, "CH2r", 6, 2},   {9, "NR", 7, 0}};
    for (auto &in : atom_dat) ff.NewAtomType(in.first, in.second, PT[in.third]);

    // Add bond types
    std::vector<stdx::quad<int, double, double, double>> bnd_dat = {
        {34, 0.198, 50181.12, 640000.0},
        {48, 0.290283, 502214.75, 2980000.0},
        {22, 0.148, 251019.84, 5730000.0},
        {49, 0.279388, 373115.59, 2390000.0},
        {50, 0.291189, 371384.73, 2190000.0},
        {10, 0.133, 417460.4, 11800000.0},
        {23, 0.148, 334693.12, 7640000.0},
        {4, 0.112, 928256.0, 37000000.0},
        {14, 0.138, 418968.0, 11000000.0},
        {31, 0.178, 376405.92, 5940000.0},
        {27, 0.153, 334748.7, 7150000.0},
        {35, 0.2, 50240.0, 628000.0},
        {6, 0.125, 418750.0, 13400000.0},
        {5, 0.123, 502282.8, 16600000.0},
        {7, 0.132, 418176.0, 12000000.0},
        {12, 0.134, 420170.4, 11700000.0},
        {39, 0.11, 292820.0, 12100000.0},
        {19, 0.143, 376670.58, 9210000.0},
        {28, 0.161, 250915.28, 4840000.0},
        {21, 0.147, 376428.78, 8710000.0},
        {37, 0.221, 52748.28, 540000.0},
        {46, 0.163299, 464531.53, 8710000.0},
        {45, 0.135, 375435.0, 10300000.0},
        {16, 0.139, 417333.6, 10800000.0},
        {52, 0.287407, 502224.92, 3040000.0},
        {41, 0.153, 376416.72, 8040000.0},
        {13, 0.136, 377318.4, 10200000.0},
        {24, 0.148, 376748.8, 8600000.0},
        {25, 0.15, 376650.0, 8370000.0},
        {17, 0.14, 334768.0, 8540000.0},
        {2, 0.1, 374000.0, 18700000.0},
        {33, 0.187, 251077.42, 3590000.0},
        {36, 0.204, 418656.96, 5030000.0},
        {15, 0.139, 334639.72, 8660000.0},
        {20, 0.1435, 251225.45, 6100000.0},
        {43, 0.176, 501811.2, 8100000.0},
        {40, 0.1758, 501907.59, 8120000.0},
        {32, 0.183, 376416.36, 5620000.0},
        {1, 0.1, 314000.0, 15700000.0},
        {3, 0.109, 292272.6, 12300000.0},
        {30, 0.178, 172360.96, 2720000.0},
        {11, 0.134, 377076.0, 10500000.0},
        {9, 0.133, 375006.8, 10600000.0},
        {47, 0.233839, 293088.43, 2680000.0},
        {26, 0.152, 250909.44, 5430000.0},
        {8, 0.133, 313802.86, 8870000.0},
        {42, 0.193799, 371824.72, 4950000.0},
        {18, 0.143, 334545.64, 8180000.0},
        {29, 0.163, 250811.36, 4720000.0},
        {38, 0.1, 464000.0, 23200000.0},
        {44, 0.1265, 419258.95, 13100000.0},
        {51, 0.2077, 342525.96, 3970000.0}};
    for (auto &dat : bnd_dat)
      ff.LinkBondTypes(
          ff.NewBondType(BondType::Harmonic, dat.first, dat.third, dat.second),
          ff.NewBondType(BondType::Quartic, dat.first, dat.fourth, dat.second));

    // Add angle types
    std::vector<stdx::quad<int, double, double, double>> ang_dat = {
        {1, 90.0, 0.11550101, 380.0},    {3, 96.0, 0.12177061, 405.0},
        {14, 109.6, 0.12142334, 450.0},  {9, 109.5, 0.08638614, 320.0},
        {40, 155.0, 0.12112698, 2215.0}, {15, 111.0, 0.14048747, 530.0},
        {29, 120.0, 0.17801113, 780.0},  {34, 125.0, 0.076490216, 375.0},
        {24, 120.0, 0.10147593, 445.0},  {37, 126.0, 0.12744672, 640.0},
        {45, 97.4, 0.14024534, 469.0},   {28, 120.0, 0.15288018, 670.0},
        {20, 116.0, 0.11421859, 465.0},  {12, 109.5, 0.12157397, 450.0},
        {27, 120.0, 0.12774922, 560.0},  {36, 126.0, 0.11448735, 575.0},
        {43, 107.57, 0.13376523, 484.0}, {46, 106.75, 0.14026005, 503.0},
        {51, 110.3, 0.14017954, 524.0},  {5, 103.0, 0.12122177, 420.0},
        {8, 109.5, 0.076912479, 285.0},  {47, 108.53, 0.12108416, 443.0},
        {19, 115.0, 0.1524165, 610.0},   {41, 180.0, 0.072640156, 91350.0},
        {16, 113.0, 0.14045138, 545.0},  {4, 100.0, 0.14008261, 475.0},
        {2, 90.0, 0.12768574, 420.0},    {50, 109.5, 0.12103261, 448.0},
        {18, 115.0, 0.11488482, 460.0},  {17, 115.0, 0.01229657, 50.0},
        {32, 123.0, 0.088743846, 415.0}, {44, 111.3, 0.16689058, 632.0},
        {39, 132.0, 0.12775497, 760.0},  {13, 109.5, 0.14052124, 520.0},
        {31, 122.0, 0.15317431, 700.0},  {53, 117.2, 0.15305438, 636.0},
        {42, 109.5, 0.11724316, 434.0},  {21, 116.0, 0.15236094, 620.0},
        {30, 121.0, 0.15312732, 685.0},  {35, 125.0, 0.1531408, 750.0},
        {6, 104.0, 0.14028506, 490.0},   {48, 109.5, 0.1670474, 618.0},
        {52, 111.4, 0.14025677, 532.0},  {54, 121.4, 0.1529482, 690.0},
        {7, 108.0, 0.12788754, 465.0},   {22, 117.0, 0.15336019, 635.0},
        {11, 109.5, 0.11480708, 425.0},  {26, 120.0, 0.12089532, 530.0},
        {10, 109.5, 0.10262668, 380.0},  {33, 124.0, 0.15266919, 730.0},
        {25, 120.0, 0.11518373, 505.0},  {38, 126.0, 0.15336544, 770.0},
        {49, 107.6, 0.14008648, 507.0},  {23, 120.0, 0.088910434, 390.0}};
    for (auto &dat : ang_dat)
      ff.LinkAngleTypes(ff.NewAngleType(AngleType::Harmonic, dat.first,
                                        dat.third, dat.second),
                        ff.NewAngleType(AngleType::CosineHarmonic, dat.first,
                                        dat.fourth, dat.second));

    // Add improper types
    std::vector<stdx::triple<int, double, double>> imp_dat = {
        {4, 0.051, 180.0},
        {2, 0.102, 35.26439},
        {5, 0.102, -35.26439},
        {3, 0.204, 0.0},
        {1, 0.051, 0.0}};
    for (auto &dat : imp_dat)
      ff.NewDihedralType(DihedralType::Improper, dat.first, dat.second,
                         dat.third);

    // Add proper types
    std::vector<stdx::quad<int, double, double, int>> prp_dat = {
        {19, 3.14, 0.0, 2},   {42, 3.5, 180.0, 2},  {17, 0.418, 0.0, 2},
        {2, 3.41, 180.0, 1},  {7, 2.79, 0.0, 1},    {24, 1.3, 0.0, 3},
        {37, 9.5, 0.0, 3},    {45, 0.4, 0.0, 6},    {15, 41.8, 180.0, 2},
        {16, 0.0, 0.0, 2},    {38, 0.0, 0.0, 4},    {41, 3.77, 0.0, 6},
        {44, 0.7, 180.0, 6},  {22, 1.05, 0.0, 3},   {6, 9.45, 180.0, 1},
        {10, 5.86, 180.0, 2}, {28, 3.65, 0.0, 3},   {13, 24.0, 180.0, 2},
        {1, 2.67, 180.0, 1},  {35, 7.69, 0.0, 3},   {11, 7.11, 180.0, 2},
        {30, 3.9, 0.0, 3},    {3, 4.97, 180.0, 1},  {20, 5.09, 0.0, 2},
        {23, 1.26, 0.0, 3},   {32, 4.69, 0.0, 3},   {26, 2.93, 0.0, 3},
        {25, 2.53, 0.0, 3},   {9, 1.53, 180.0, 2},  {12, 16.7, 180.0, 2},
        {27, 3.19, 0.0, 3},   {14, 33.5, 180.0, 2}, {8, 5.35, 0.0, 1},
        {33, 5.44, 0.0, 3},   {34, 5.92, 0.0, 3},   {5, 9.35, 180.0, 1},
        {29, 3.77, 0.0, 3},   {39, 1.0, 180.0, 6},  {4, 5.86, 180.0, 1},
        {31, 4.18, 0.0, 3},   {43, 2.8, 0.0, 3},    {36, 8.62, 0.0, 3},
        {40, 1.0, 0.0, 6},    {21, 16.7, 0.0, 2},   {18, 2.09, 0.0, 2}};
    for (auto &dat : prp_dat)
      ff.NewDihedralType(DihedralType::Proper, dat.first, dat.second, dat.third,
                         dat.fourth);

    return ff;
  }
} // namespace indigox
