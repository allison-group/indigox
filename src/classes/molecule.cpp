#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/serialise.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid molecule instance")
#else
#define _sanity_check_(x)
#endif

namespace indigox {

  using Data = CalculatedData;

  // =======================================================================
  // == SERIALISATION ======================================================
  // =======================================================================

  template <typename Archive>
  void Molecule::Impl::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("name", name),
            INDIGOX_SERIAL_NVP("next_uid", next_unique_id),
            INDIGOX_SERIAL_NVP("charge", molecular_charge),
            INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("bonds", bonds),
            INDIGOX_SERIAL_NVP("angles", angles),
            INDIGOX_SERIAL_NVP("dihedrals", dihedrals),
            INDIGOX_SERIAL_NVP("forcefield", forcefield),
            INDIGOX_SERIAL_NVP("graph", molecular_graph),
            INDIGOX_SERIAL_NVP("formula", cached_formula),
            INDIGOX_SERIAL_NVP("calculated", calculated_data));
  }

  template <typename Archive>
  void Molecule::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(Molecule);

  // =======================================================================
  // == OPERATORS ==========================================================
  // =======================================================================

  bool Molecule::operator==(const Molecule &mol) const {
    _sanity_check_(*this);
    _sanity_check_(mol);
    return m_data == mol.m_data;
  }

  bool Molecule::operator<(const Molecule &mol) const {
    _sanity_check_(*this);
    _sanity_check_(mol);
    return m_data < mol.m_data;
  }

  bool Molecule::operator>(const Molecule &mol) const {
    _sanity_check_(*this);
    _sanity_check_(mol);
    return m_data > mol.m_data;
  }

  // =======================================================================
  // == CONSTRUCTION =======================================================
  // =======================================================================

  Molecule::Impl::Impl(std::string n)
      : name(n), next_unique_id(0), molecular_charge(0), calculated_data(0),
        cached_formula("") {}

  Molecule::Molecule(std::string n) : m_data(std::make_shared<Impl>(n)) {
    m_data->molecular_graph = graph::MolecularGraph(*this);
  }

  // =======================================================================
  // == STATE CHECKING =====================================================
  // =======================================================================

  int64_t Molecule::Impl::FindBond(const Atom &a, const Atom &b) const {
    //    _sanity_check_(*this);
    int64_t pos = 0;
    for (const Bond &bnd : a.GetBonds()) {
      if (bnd.GetAtoms()[0] == b || bnd.GetAtoms()[1] == b) break;
      ++pos;
    }
    return pos == a.NumBonds() ? -1 : pos;
  }

  int64_t Molecule::Impl::FindAngle(const Atom &a, const Atom &b,
                                    const Atom &c) const {
    //    _sanity_check_(*this);
    int64_t pos = 0;
    for (const Angle &ang : b.GetAngles()) {
      if ((ang.GetAtoms()[0] == a && ang.GetAtoms()[2] == c) ||
          (ang.GetAtoms()[2] == a && ang.GetAtoms()[1] == c)) {
        break;
      }
      ++pos;
    }
    return pos == b.NumAngles() ? -1 : pos;
  }

  int64_t Molecule::Impl::FindDihedral(const Atom &a, const Atom &b,
                                       const Atom &c, const Atom &d) const {
    //    _sanity_check_(*this);
    int64_t pos = 0;
    for (const Dihedral &dhd : a.GetDihedrals()) {
      if ((dhd.GetAtoms()[0] == a && dhd.GetAtoms()[1] == b &&
           dhd.GetAtoms()[2] == c && dhd.GetAtoms()[3] == d) ||
          (dhd.GetAtoms()[3] == a && dhd.GetAtoms()[2] == b &&
           dhd.GetAtoms()[1] == c && dhd.GetAtoms()[0] == d)) {
        break;
      }
      ++pos;
    }
    return pos == a.NumDihedrals() ? -1 : pos;
  }

  bool Molecule::HasAtom(const Atom &atom) const {
    _sanity_check_(*this);
    return atom.GetMolecule() == *this;
  }

  bool Molecule::HasBond(const Bond &bond) const {
    _sanity_check_(*this);
    return bond.GetMolecule() == *this;
  }

  bool Molecule::HasBond(const Atom &a, const Atom &b) const {
    _sanity_check_(*this);
    return (HasAtom(a) && HasAtom(b)) ? (m_data->FindBond(a, b) != -1) : false;
  }

  bool Molecule::HasAngle(const Angle &angle) const {
    _sanity_check_(*this);
    return angle.GetMolecule() == *this;
  }

  bool Molecule::HasAngle(const Atom &a, const Atom &b, const Atom &c) {
    _sanity_check_(*this);
    PerceiveAngles();
    return (HasAtom(a) && HasAtom(b) && HasAtom(c))
               ? (m_data->FindAngle(a, b, c) != -1)
               : false;
  }

  bool Molecule::HasDihedral(const Dihedral &dihedral) const {
    _sanity_check_(*this);
    return dihedral.GetMolecule() == *this;
  }

  bool Molecule::HasDihedral(const Atom &a, const Atom &b, const Atom &c,
                             const Atom &d) {
    _sanity_check_(*this);
    PerceiveDihedrals();
    return (HasAtom(a) && HasAtom(b) && HasAtom(c) && HasAtom(d))
               ? (m_data->FindDihedral(a, b, c, d) != -1)
               : false;
  }

  bool Molecule::HasForcefield() const {
    _sanity_check_(*this);
    return bool(m_data->forcefield);
  }

  int64_t Molecule::NumAtoms() const {
    _sanity_check_(*this);
    return static_cast<int64_t>(m_data->atoms.size());
  }

  int64_t Molecule::NumBonds() const {
    _sanity_check_(*this);
    return static_cast<int64_t>(m_data->bonds.size());
  }

  int64_t Molecule::NumAngles() {
    _sanity_check_(*this);
    PerceiveAngles();
    return static_cast<int64_t>(m_data->angles.size());
  }

  int64_t Molecule::NumDihedrals() {
    _sanity_check_(*this);
    PerceiveDihedrals();
    return static_cast<int64_t>(m_data->dihedrals.size());
  }

  void Molecule::ReserveAtoms(int64_t num) {
    _sanity_check_(*this);
    if (num > 0) { m_data->atoms.reserve(num); }
  }

  void Molecule::ReserveBonds(int64_t num) {
    _sanity_check_(*this);
    if (num > 0) { m_data->bonds.reserve(num); }
  }

  // =======================================================================
  // == STATE GETTING ======================================================
  // =======================================================================

  Atom Molecule::GetAtom(uint32_t pos) const {
    _sanity_check_(*this);
    return (pos < NumAtoms()) ? m_data->atoms[pos] : Atom();
  }

  Atom Molecule::GetAtomID(int64_t id) const {
    _sanity_check_(*this);
    for (const Atom &atm : m_data->atoms) {
      if (atm.GetID() == id) { return atm; }
    }
    return Atom();
  }

  Atom Molecule::GetAtomTag(int64_t tag) const {
    _sanity_check_(*this);
    for (const Atom &atm : m_data->atoms) {
      if (atm.GetTag() == tag) { return atm; }
    }
    return Atom();
  }

  Bond Molecule::GetBond(uint32_t pos) const {
    _sanity_check_(*this);
    return (pos < NumBonds()) ? m_data->bonds[pos] : Bond();
  }

  Bond Molecule::GetBond(const Atom &a, const Atom &b) const {
    _sanity_check_(*this);
    int64_t pos = m_data->FindBond(a, b);
    return (HasAtom(a) && pos != -1) ? a.GetBonds()[pos] : Bond();
  }

  Bond Molecule::GetBondID(int64_t id) const {
    _sanity_check_(*this);
    for (const Bond &bnd : m_data->bonds) {
      if (bnd.GetID() == id) { return bnd; }
    }
    return Bond();
  }

  Bond Molecule::GetBondTag(int64_t tag) const {
    _sanity_check_(*this);
    for (const Bond &bnd : m_data->bonds) {
      if (bnd.GetTag() == tag) { return bnd; }
    }
    return Bond();
  }

  Angle Molecule::GetAngle(uint32_t pos) {
    _sanity_check_(*this);
    PerceiveAngles();
    return (pos < NumAngles()) ? m_data->angles[pos] : Angle();
  }

  Angle Molecule::GetAngle(const Atom &a, const Atom &b, const Atom &c) {
    _sanity_check_(*this);
    PerceiveAngles();
    int64_t pos = m_data->FindAngle(a, b, c);
    return (HasAtom(b) && pos != -1) ? b.GetAngles()[pos] : Angle();
  }

  Angle Molecule::GetAngleID(int64_t id) const {
    _sanity_check_(*this);
    for (const Angle &ang : m_data->angles) {
      if (ang.GetID() == id) { return ang; }
    }
    return Angle();
  }

  Angle Molecule::GetAngleTag(int64_t tag) const {
    _sanity_check_(*this);
    for (const Angle &ang : m_data->angles) {
      if (ang.GetTag() == tag) { return ang; }
    }
    return Angle();
  }

  Dihedral Molecule::GetDihedral(uint32_t pos) {
    _sanity_check_(*this);
    PerceiveDihedrals();
    return (pos < NumDihedrals()) ? m_data->dihedrals[pos] : Dihedral();
  }

  Dihedral Molecule::GetDihedral(const Atom &a, const Atom &b, const Atom &c,
                                 const Atom &d) {
    _sanity_check_(*this);
    PerceiveDihedrals();
    int64_t pos = m_data->FindDihedral(a, b, c, d);
    return (HasAtom(a) && pos != -1) ? a.GetDihedrals()[pos] : Dihedral();
  }

  Dihedral Molecule::GetDihedralID(int64_t id) const {
    _sanity_check_(*this);
    for (const Dihedral &ang : m_data->dihedrals) {
      if (ang.GetID() == id) { return ang; }
    }
    return Dihedral();
  }

  Dihedral Molecule::GetDihedralTag(int64_t tag) const {
    _sanity_check_(*this);
    for (const Dihedral &ang : m_data->dihedrals) {
      if (ang.GetTag() == tag) { return ang; }
    }
    return Dihedral();
  }

  Residue Molecule::GetResidueID(int32_t id) { return GetResidues()[id]; }

  std::string Molecule::GetFormula() {
    _sanity_check_(*this);
    if (!m_data->Test(Data::Formula)) {
      std::map<std::string, size_t> e_count;
      for (const Atom &atm : m_data->atoms)
        e_count[atm.GetElement().GetSymbol()]++;
      std::stringstream ss;
      if (e_count["C"]) ss << "C";
      if (e_count["C"] > 1) ss << e_count["C"];
      if (e_count["H"]) ss << "H";
      if (e_count["H"] > 1) ss << e_count["H"];
      for (auto &e : e_count) {
        if (e.first != "C" && e.first != "H") {
          ss << e.first;
          if (e.second > 1) ss << e.second;
        }
      }
      m_data->Set(Data::Formula);
      m_data->cached_formula = ss.str();
    }
    return m_data->cached_formula;
  }

  const graph::MolecularGraph &Molecule::GetGraph() const {
    _sanity_check_(*this);
    return m_data->molecular_graph;
  }

  const graph::CondensedMolecularGraph &Molecule::GetCondensedGraph() const {
    _sanity_check_(*this);
    if (!m_data->Test(Data::CondensedGraph)) {
      m_data->condensed_molecular_graph =
          graph::Condense(m_data->molecular_graph);
      m_data->Set(Data::CondensedGraph);
    }
    return m_data->condensed_molecular_graph;
  }

  const std::string &Molecule::GetName() const {
    _sanity_check_(*this);
    return m_data->name;
  }

  int32_t Molecule::GetMolecularCharge() const {
    _sanity_check_(*this);
    return m_data->molecular_charge;
  }

  const Molecule::MoleculeAtoms &Molecule::GetAtoms() const {
    _sanity_check_(*this);
    return m_data->atoms;
  }

  const Molecule::MoleculeBonds &Molecule::GetBonds() const {
    _sanity_check_(*this);
    return m_data->bonds;
  }

  const Molecule::MoleculeAngles &Molecule::GetAngles() {
    _sanity_check_(*this);
    PerceiveAngles();
    return m_data->angles;
  }

  const Molecule::MoleculeDihedrals &Molecule::GetDihedrals() {
    _sanity_check_(*this);
    PerceiveDihedrals();
    return m_data->dihedrals;
  }

  const Molecule::MoleculeResidues &Molecule::GetResidues() {
    _sanity_check_(*this);
    PerceiveResidues();
    return m_data->residues;
  }

  const Forcefield &Molecule::GetForcefield() const {
    _sanity_check_(*this);
    return m_data->forcefield;
  }

  // =======================================================================
  // == STATE MODIFYING ====================================================
  // =======================================================================

  void Molecule::ReorderAtoms(MoleculeAtoms &new_order) {
    eastl::vector_set<Atom> all_atoms(m_data->atoms.begin(),
                                      m_data->atoms.end());
    MoleculeAtoms swap_order(new_order);
    for (Atom atm : new_order) {
      if (!HasAtom(atm)) {
        throw std::runtime_error(
            "New atom order contains atoms not part of this molecule");
      }
      all_atoms.erase(atm);
    }
    swap_order.insert(swap_order.end(), all_atoms.begin(), all_atoms.end());
    m_data->atoms = swap_order;
  }

  void Molecule::OptimiseChargeGroups() {
    std::vector<std::vector<Atom>> charge_groups =
        algorithm::OptimalChargeGroups(*this);
    MoleculeAtoms new_order;
    int32_t charge_group = -1;
    for (std::vector<Atom> grp : charge_groups) {
      ++charge_group;
      new_order.insert(new_order.end(), grp.begin(), grp.end());
      for (Atom atm : grp) atm.SetChargeGroupID(charge_group);
    }
    ReorderAtoms(new_order);
  }

  Atom Molecule::NewAtom() {
    return NewAtom(GetPeriodicTable().GetUndefined(), 0.0, 0.0, 0.0);
  }

  Atom Molecule::NewAtom(const Element &element) {
    return NewAtom(element, 0.0, 0.0, 0.0);
  }

  Atom Molecule::NewAtom(const Element &element, double x, double y, double z) {
    _sanity_check_(*this);
    m_data->Reset();
    Atom atom = Atom(*this, element, x, y, z, "");
    atom.m_data->unique_id = m_data->next_unique_id++;
    m_data->atoms.emplace_back(atom);
    m_data->molecular_graph.AddVertex(atom);
    return atom;
  }

  Bond Molecule::NewBond(const Atom &a, const Atom &b) {
    _sanity_check_(*this);
    Bond bnd;

    if (HasAtom(a) && HasAtom(b)) {
      int64_t pos = m_data->FindBond(a, b);
      if (pos != -1) {
        bnd = a.GetBonds()[pos];
      } else {
        m_data->Reset();
        bnd = Bond(a, b, *this, BondOrder::SINGLE);
        bnd.m_data->unique_id = m_data->next_unique_id++;
        bnd.m_data->atoms[0].AddBond(bnd);
        bnd.m_data->atoms[1].AddBond(bnd);
        m_data->molecular_graph.AddEdge(bnd);
        m_data->bonds.emplace_back(bnd);
      }
    }

    return bnd;
  }

  Angle Molecule::NewAngle(const Atom &a, const Atom &b, const Atom &c) {
    // No need for logic checks as no ability for user to add angles
    _sanity_check_(*this);
    Angle ang = Angle(a, b, c, *this);
    ang.m_data->unique_id = m_data->next_unique_id++;
    ang.m_data->atoms[0].AddAngle(ang);
    ang.m_data->atoms[1].AddAngle(ang);
    ang.m_data->atoms[2].AddAngle(ang);
    m_data->angles.emplace_back(ang);
    return ang;
  }

  Dihedral Molecule::NewDihedral(const Atom &a, const Atom &b, const Atom &c,
                                 const Atom &d, bool manual) {
    Dihedral dhd;
    if (manual) {
      // Run sanity checks on manual
      _sanity_check_(*this);
      if (a == b || a == c || a == d || b == c || b == d || c == d) {
        return dhd;
      }
      if (!(HasAtom(a) && HasAtom(b) && HasAtom(c) && HasAtom(d))) {
        return dhd;
      }
    }

    int64_t pos = m_data->FindDihedral(a, b, c, d);
    if (pos != -1) {
      dhd = a.GetDihedrals()[pos];
    } else {
      dhd = Dihedral(a, b, c, d, *this);
      dhd.m_data->unique_id = m_data->next_unique_id++;
      dhd.m_data->atoms[0].AddDihedral(dhd);
      dhd.m_data->atoms[1].AddDihedral(dhd);
      dhd.m_data->atoms[2].AddDihedral(dhd);
      dhd.m_data->atoms[3].AddDihedral(dhd);
      m_data->dihedrals.emplace_back(dhd);
    }
    return dhd;
  }

  Dihedral Molecule::NewDihedral(const Atom &a, const Atom &b, const Atom &c,
                                 const Atom &d) {
    return NewDihedral(a, b, c, d, true);
  }

  bool Molecule::RemoveAtom(const Atom &atom) {
    _sanity_check_(*this);
    if (!HasAtom(atom)) return false;
    m_data->Reset();
    // Remove all bonds this atom is part of from molecule
    auto bnd_pred = [&atom](Bond bnd) { // Predicate checks if atom in bnd
      return !(bnd.m_data->atoms[0] == atom || bnd.m_data->atoms[1] == atom);
    };
    // Partition so can remove bonds from other atoms
    auto bnd_pos =
        std::partition(m_data->bonds.begin(), m_data->bonds.end(), bnd_pred);
    auto bnd_pos_erase = bnd_pos;
    for (auto it = m_data->bonds.end(); bnd_pos != it; ++bnd_pos) {
      Bond bnd = *bnd_pos;
      bnd.m_data->atoms[0].RemoveBond(bnd);
      bnd.m_data->atoms[1].RemoveBond(bnd);
      bnd.Reset();
    }
    m_data->bonds.erase(bnd_pos_erase, m_data->bonds.end());

    // Remove all angles this atom is part of from molecule
    auto ang_pred = [&atom](Angle ang) {
      return !(ang.m_data->atoms[0] == atom || ang.m_data->atoms[1] == atom ||
               ang.m_data->atoms[2] == atom);
    };
    auto ang_pos =
        std::partition(m_data->angles.begin(), m_data->angles.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = m_data->angles.end(); ang_pos != it; ++ang_pos) {
      Angle ang = *ang_pos;
      ang.m_data->atoms[0].RemoveAngle(ang);
      ang.m_data->atoms[1].RemoveAngle(ang);
      ang.m_data->atoms[2].RemoveAngle(ang);
      ang.Reset();
    }
    m_data->angles.erase(ang_pos_erase, m_data->angles.end());

    // Remove all dihedrals this atom is part of from molecule
    auto dhd_pred = [&atom](Dihedral dhd) {
      return !(dhd.m_data->atoms[0] == atom || dhd.m_data->atoms[1] == atom ||
               dhd.m_data->atoms[2] == atom || dhd.m_data->atoms[3] == atom);
    };
    auto dhd_pos = std::partition(m_data->dihedrals.begin(),
                                  m_data->dihedrals.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = m_data->dihedrals.end(); dhd_pos != it; ++dhd_pos) {
      Dihedral dhd = *dhd_pos;
      dhd.m_data->atoms[0].RemoveDihedral(dhd);
      dhd.m_data->atoms[1].RemoveDihedral(dhd);
      dhd.m_data->atoms[2].RemoveDihedral(dhd);
      dhd.m_data->atoms[3].RemoveDihedral(dhd);
      dhd.Reset();
    }
    m_data->dihedrals.erase(dhd_pos_erase, m_data->dihedrals.end());

    // Remove the atom from the molecule
    Atom atm = atom;
    graph::MGVertex v = m_data->molecular_graph.GetVertex(atm);
    m_data->molecular_graph.RemoveVertex(v);
    m_data->atoms.erase(
        std::find(m_data->atoms.begin(), m_data->atoms.end(), atm));
    atm.Reset();
    return true;
  }

#define Square(val) val *val

  bool Molecule::RemoveBond(const Bond &bond) {
    _sanity_check_(*this);
    if (!HasBond(bond)) return false;
    m_data->Reset();
    Atom a = bond.GetAtoms()[0];
    Atom b = bond.GetAtoms()[1];
    // Remove the bond from each atom
    a.RemoveBond(bond);
    b.RemoveBond(bond);

    // Remove all angles this bond is part of from molecule
    auto ang_pred = [&](Angle ang) {
      int64_t b1_pos =
          m_data->FindBond(ang.m_data->atoms[1], ang.m_data->atoms[0]);
      int64_t b2_pos =
          m_data->FindBond(ang.m_data->atoms[1], ang.m_data->atoms[2]);
      return (ang.GetAtoms()[1].GetBonds()[b1_pos] != bond &&
              ang.GetAtoms()[1].GetBonds()[b2_pos] != bond);
    };
    auto ang_pos =
        std::partition(m_data->angles.begin(), m_data->angles.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = m_data->angles.end(); it != ang_pos; ++ang_pos) {
      Angle ang = *ang_pos;
      ang.m_data->atoms[0].RemoveAngle(ang);
      ang.m_data->atoms[1].RemoveAngle(ang);
      ang.m_data->atoms[2].RemoveAngle(ang);
      ang.Reset();
    }
    m_data->angles.erase(ang_pos_erase, m_data->angles.end());

    // Remove all dihedrals this bond is part of from molecule
    auto dhd_pred = [&](Dihedral dhd) {
      int64_t b1_pos = m_data->FindBond(dhd.GetAtoms()[0], dhd.GetAtoms()[1]);
      int64_t b2_pos = m_data->FindBond(dhd.GetAtoms()[1], dhd.GetAtoms()[2]);
      int64_t b3_pos = m_data->FindBond(dhd.GetAtoms()[2], dhd.GetAtoms()[3]);
      if (b1_pos == -1 || b2_pos == -1 || b3_pos == -1) { return true; }
      Atom a = dhd.GetAtoms()[0];
      return (a.GetBonds()[b1_pos] != bond && a.GetBonds()[b2_pos] != bond &&
              a.GetBonds()[b3_pos] != bond);
    };
    auto dhd_pos = std::partition(m_data->dihedrals.begin(),
                                  m_data->dihedrals.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = m_data->dihedrals.end(); dhd_pos != it; ++dhd_pos) {
      Dihedral dhd = *dhd_pos;
      dhd.m_data->atoms[0].RemoveDihedral(dhd);
      dhd.m_data->atoms[1].RemoveDihedral(dhd);
      dhd.m_data->atoms[2].RemoveDihedral(dhd);
      dhd.m_data->atoms[3].RemoveDihedral(dhd);
      dhd.Reset();
    }
    m_data->dihedrals.erase(dhd_pos_erase, m_data->dihedrals.end());

    // Remove the bond from the molecule
    Bond bnd = bond;
    graph::MGEdge e = m_data->molecular_graph.GetEdge(bnd);
    m_data->molecular_graph.RemoveEdge(e);
    m_data->bonds.erase(
        std::find(m_data->bonds.begin(), m_data->bonds.end(), bond));
    bnd.Reset();
    return true;
  }

  bool Molecule::RemoveBond(const Atom &a, const Atom &b) {
    return RemoveBond(GetBond(a, b));
  }

  int64_t Molecule::PerceiveAngles() {
    _sanity_check_(*this);
    if (m_data->Test(Data::AnglePerception)) { return 0; }
    m_data->Set(Data::AnglePerception);

    // Expected number of angles
    auto sum = [&](size_t current, Atom v) -> size_t {
      size_t degree = v.NumBonds();
      if (degree < 2) return current;
      return current + degree * (degree - 1) / 2;
    };
    size_t count =
        std::accumulate(m_data->atoms.begin(), m_data->atoms.end(), 0, sum);
    m_data->angles.reserve(count);
    count = 0;

    std::vector<Atom> nbrs;
    nbrs.reserve(10);
    // Adding new angles
    for (Atom at : m_data->atoms) {
      if (at.NumBonds() < 2) { continue; }
      for (Bond bn : at.m_data->bonds) {
        if (bn.GetAtoms()[0] == at) {
          nbrs.emplace_back(bn.GetAtoms()[1]);
        } else {
          nbrs.emplace_back(bn.GetAtoms()[0]);
        }
      }

      for (size_t i = 0; i < nbrs.size() - 1; ++i) {
        for (size_t j = i + 1; j < nbrs.size(); ++j) {
          if (m_data->FindAngle(nbrs[i], at, nbrs[j]) != -1) { continue; }
          NewAngle(nbrs[i], at, nbrs[j]);
          ++count;
        }
      }
      nbrs.clear();
    }
    return count;
  }

  int64_t Molecule::PerceiveDihedrals() {
    if (m_data->Test(Data::DihedralPerception)) { return 0; }
    m_data->Set(Data::DihedralPerception);

    // Expected number of dihedrals
    auto sum = [&](size_t current, Bond b) -> size_t {
      size_t b_degree = b.GetAtoms()[0].NumBonds();
      if (b_degree < 2) { return current; }
      size_t c_degree = b.GetAtoms()[1].NumBonds();
      if (c_degree < 2) { return current; }
      return current + (b_degree - 1) * (c_degree - 1);
    };
    size_t count =
        std::accumulate(m_data->bonds.begin(), m_data->bonds.end(), 0, sum);
    m_data->dihedrals.reserve(count);
    count = 0;

    // Adding new dihedrals
    std::vector<Atom> B_nbrs, C_nbrs;
    B_nbrs.reserve(10);
    C_nbrs.reserve(10);
    for (Bond bn : m_data->bonds) {
      Atom B = bn.GetAtoms()[0];
      Atom C = bn.GetAtoms()[1];
      if (B.NumBonds() < 2 || C.NumBonds() < 2) { continue; }

      for (Bond bn : B.GetBonds()) {
        if (bn.GetAtoms()[0] == B) {
          B_nbrs.emplace_back(bn.GetAtoms()[1]);
        } else {
          B_nbrs.emplace_back(bn.GetAtoms()[0]);
        }
      }
      for (Bond bn : C.GetBonds()) {
        if (bn.GetAtoms()[0] == C) {
          C_nbrs.emplace_back(bn.GetAtoms()[1]);
        } else {
          C_nbrs.emplace_back(bn.GetAtoms()[0]);
        }
      }

      for (size_t i = 0; i < B_nbrs.size(); ++i) {
        if (B_nbrs[i] == C) { continue; }
        for (size_t j = 0; j < C_nbrs.size(); ++j) {
          if (C_nbrs[j] == B) { continue; }
          if (m_data->FindDihedral(B_nbrs[i], B, C, C_nbrs[j]) != -1) continue;
          NewDihedral(B_nbrs[i], B, C, C_nbrs[j], false);
          ++count;
        }
      }
      B_nbrs.clear();
      C_nbrs.clear();
    }

    // Determining priorities
    for (Bond bn : m_data->bonds) {
      Atom B = bn.GetAtoms()[0];
      Atom C = bn.GetAtoms()[1];
      if (B.NumBonds() < 2 || C.NumBonds() < 2) { continue; }
      std::vector<Dihedral> bnd_dhds;
      for (Dihedral dhd : B.GetDihedrals()) {
        Atom b = dhd.GetAtoms()[1];
        Atom c = dhd.GetAtoms()[2];
        if ((b == B && c == C) || (b == C && c == B)) {
          bnd_dhds.emplace_back(dhd);
        }
      }
      std::sort(bnd_dhds.begin(), bnd_dhds.end(), [](Dihedral a, Dihedral b) {
        Atom a1 = a.GetAtoms()[0];
        Atom a4 = a.GetAtoms()[3];
        Atom b1 = b.GetAtoms()[0];
        Atom b4 = b.GetAtoms()[3];
        int32_t w_a = Square(a1.GetElement().GetAtomicNumber()) +
                      Square(a4.GetElement().GetAtomicNumber());
        int32_t w_b = Square(b1.GetElement().GetAtomicNumber()) +
                      Square(b4.GetElement().GetAtomicNumber());
        return w_a < w_b;
      });
      for (uint32_t i = 0; i < bnd_dhds.size(); ++i) {
        // Only worry about the H's for now.
        if (bnd_dhds[i].GetAtoms()[0].GetElement() != "H" &&
            bnd_dhds[i].GetAtoms()[3].GetElement() != "H") {
          continue;
        }
        bnd_dhds[i].m_data->priority = i;
      }
    }
    return count;
  }

  // Finds the edges of a graph corresponding to peptide bonds in a molecule
  void PeptideBonds(const graph::MolecularGraph &G,
                    std::vector<graph::MGEdge> &edges) {
    for (graph::MGEdge e : G.GetEdges()) {
      if (!e.GetBond().IsAmideBond()) continue;
      for (Atom atm : e.GetBond().GetAtoms()) {
        if (atm.GetElement() == "N" &&
            atm.NumHydrogenBonds() + atm.GetImplicitCount() <= 1) {
          edges.push_back(e);
        }
      }
    }
  }

  int32_t Molecule::PerceiveResidues() {
    using namespace graph;
    _sanity_check_(*this);
    if (m_data->Test(Data::ResiduePerception)) {
      return (int32_t)m_data->residues.size();
    }
    m_data->Set(Data::ResiduePerception);

    MolecularGraph graph = m_data->molecular_graph;

    MolecularGraph::VertContain all_vertices = graph.GetVertices();
    MolecularGraph residue_graph = graph.Subgraph(all_vertices);

    // Identify peptide bonds
    std::vector<graph::MGEdge> to_remove;
    PeptideBonds(residue_graph, to_remove);
    //! \todo Other types of bonds to break on for residue identification?

    // Remove the peptide bonds to get the components of the graph as residues
    std::vector<MGVertex> res_vert;
    res_vert.reserve(2 * to_remove.size());
    for (MGEdge e : to_remove) {
      residue_graph.RemoveEdge(e);
      res_vert.emplace_back(graph.GetSourceVertex(e));
      res_vert.emplace_back(graph.GetTargetVertex(e));
    }
    // Uniquify the residue vertices, just in case
    std::sort(res_vert.begin(), res_vert.end());
    auto last = std::unique(res_vert.begin(), res_vert.end());
    res_vert.erase(last, res_vert.end());

    // Make a graph of only the vertices of residue breaking bonds
    MolecularGraph residue_ordering = graph.Subgraph(res_vert, to_remove);
    // Add bonds between every atom pair within a component
    for (auto component : residue_graph.GetConnectedComponents()) {
      for (MGVertex source : component) {
        if (!residue_ordering.HasVertex(source)) continue;
        for (MGVertex target : component) {
          if (!residue_ordering.HasVertex(target)) continue;
          Bond tmp_bnd(source.GetAtom(), target.GetAtom(), *this,
                       BondOrder::SINGLE);
          residue_ordering.AddEdge(tmp_bnd);
        }
      }
    }

    // order vertices based on degree and element. carbons and d(1) first.
    std::sort(res_vert.begin(), res_vert.end(),
              [&residue_ordering](MGVertex u, MGVertex v) {
                if (residue_ordering.Degree(u) != residue_ordering.Degree(v))
                  return residue_ordering.Degree(u) <
                         residue_ordering.Degree(v);
                return u.GetAtom().GetElement().GetAtomicNumber() <
                       v.GetAtom().GetElement().GetAtomicNumber();
              });

    //----
    //- Order vertices to get sensible order for component storing
    //---

    // For every component, Do a depth first search starting from the first C
    std::vector<MGVertex> ordered_vertices;
    for (auto &component : residue_ordering.GetConnectedComponents()) {
      MGVertex source;
      for (MGVertex test : res_vert) {
        if (std::find(component.begin(), component.end(), test) ==
            component.end())
          continue;
        source = test;
        break;
      }
      if (!source) throw std::runtime_error("Something went wrong");

      algorithm::TraversalResults<MGVertex> dfs =
          algorithm::DepthFirstSearch(residue_ordering, source);
      std::pair<MGVertex, int32_t> l_path =
          std::make_pair(dfs.furthest, dfs.path_lengths[dfs.furthest]);

      eastl::vector_set<MGVertex> seen;
      while (component.size() > seen.size()) {
        std::vector<MGVertex> current_path;
        current_path.reserve(l_path.second);
        while (l_path.first && seen.find(l_path.first) == seen.end()) {
          current_path.push_back(l_path.first);
          seen.insert(l_path.first);
          l_path.first = dfs.predecessors[l_path.first];
        }
        ordered_vertices.insert(ordered_vertices.end(), current_path.rbegin(),
                                current_path.rend());
        l_path.second = -1;
        for (auto &vl : dfs.path_lengths) {
          if (seen.find(vl.first) != seen.end()) continue;
          if (vl.second > l_path.second)
            l_path = std::make_pair(vl.first, vl.second);
        }
      }
    }

    MolecularGraph::ComponentContain components =
        residue_graph.GetConnectedComponents();
    // Figure out the order for the components
    std::vector<int32_t> order;
    for (MGVertex v : ordered_vertices) {
      int32_t index = -1;
      for (auto &component : components) {
        ++index;
        if (std::find(order.begin(), order.end(), index) != order.end())
          continue;
        if (std::find(component.begin(), component.end(), v) !=
            component.end()) {
          order.push_back(index);
          break;
        }
      }
    }

    MoleculeAtoms new_order;
    int32_t res_id = 1;
    m_data->residues.clear();
    for (int32_t index : order) {
      MolecularGraph graph = residue_graph.Subgraph(components[index]);
      // order component vertices. dfs from first in ordered_vertices
      MGVertex source;
      eastl::vector_set<MGVertex> comp_vert(graph.GetVertices().begin(),
                                            graph.GetVertices().end());
      for (MGVertex v : ordered_vertices) {
        if (comp_vert.find(v) != comp_vert.end()) {
          source = v;
          break;
        }
      }
      if (!source) throw std::runtime_error("Something went wrong");
      auto dfs = algorithm::DepthFirstSearch(graph, source);

      std::vector<MGVertex> order;
      eastl::vector_set<MGVertex> seen;
      MGVertex start = dfs.furthest;
      while (order.size() < dfs.discover_order.size()) {
        std::vector<MGVertex> path_order;
        while (start && seen.find(start) == seen.end()) {
          for (MGVertex nbr : graph.GetNeighbours(start)) {
            if (graph.Degree(nbr) != 1) continue;
            if (seen.find(nbr) == seen.end()) {
              seen.insert(nbr);
              path_order.push_back(nbr);
            }
          }
          path_order.push_back(start);
          seen.insert(start);
          start = dfs.predecessors[start];
        }
        order.insert(order.end(), path_order.rbegin(), path_order.rend());
        int32_t furthest_length = -1;
        for (MGVertex v : dfs.discover_order) {
          if (seen.find(v) != seen.end()) continue;
          if (dfs.path_lengths[v] > furthest_length) {
            start = v;
            furthest_length = dfs.path_lengths[v];
          }
        }
      }

      std::vector<Atom> ordered_atoms;
      ordered_atoms.reserve(order.size());
      for (MGVertex v : order) {
        Atom atm = v.GetAtom();
        new_order.emplace_back(atm);
        ordered_atoms.emplace_back(atm);
        atm.m_data->residue_id = res_id;
        atm.m_data->residue_name = "RS" + std::to_string(res_id);
      }
      Residue res(ordered_atoms, *this);
      m_data->residues.emplace_back(res);

      ++res_id;
    }
    ReorderAtoms(new_order);

    return (int32_t)m_data->residues.size();
  }

  void Molecule::ModificationMade() { m_data->Reset(); }

  // =======================================================================
  // == STATE SETTING ======================================================
  // =======================================================================

  void Molecule::SetName(std::string name) {
    _sanity_check_(*this);
    m_data->name = name;
  }

  void Molecule::SetMolecularCharge(int32_t q) {
    _sanity_check_(*this);
    m_data->Reset();
    m_data->molecular_charge = q;
  }

  void Molecule::SetForcefield(const Forcefield &ff) {
    _sanity_check_(*this);
    if (HasForcefield()) {
      throw std::runtime_error("Molecule already has a forcefield set.");
    }
    m_data->forcefield = ff;
  }

  void Molecule::ResetForcefield(const Forcefield &ff) {
    _sanity_check_(*this);
    m_data->forcefield = ff;
  }

  // =======================================================================
  // == OUTPUTTING =========================================================
  // =======================================================================

  void SaveMolecule(const Molecule &mol, std::string path) {
    using Archive = cereal::PortableBinaryOutputArchive;
    std::ofstream os(path);
    if (!os.is_open()) throw std::runtime_error("Unable to open output stream");
    Archive archive(os);
    std::string stype("Molecule");
    archive(stype, mol);
  }

  Molecule LoadMolecule(std::string path) {
    using Archive = cereal::PortableBinaryInputArchive;
    std::ifstream is(path);
    if (!is.is_open()) throw std::runtime_error("Unable to open input stream");
    std::string stype;
    Archive archive(is);
    archive(stype);
    if (stype != "Molecule") throw std::runtime_error("Not a Molecule file");
    Molecule mol;
    archive(mol);
    return mol;
  }

} // namespace indigox
