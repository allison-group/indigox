#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/doctest_proxy.hpp>
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

  test_suite_open("Atom");

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
      : molecule(m), element(e), formal_charge(0), implicit_hydrogens(0),
        position(x, y, z), name(n), partial_charge(0.0),
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
    return m_data->element;
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

  /*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXAtom serialisation", T,
    ixatom_serial) { using In = typename T::t1; using Out = typename
    cereal::traits::detail::get_output_from_input<In>::type;
      test::AtomTestFixture fixture;
      Atom saved = fixture.atm.imp;
      Bond bond = test::CreateGenericTestBond().imp;
      Angle angle = test::CreateGenericTestAngle().imp;
      Dihedral dihedral = test::CreateGenericTestDihedral().imp;
      saved->SetElement("Tc");
      saved->SetTag(3);
      saved->SetImplicitCount(23);
      saved->SetName("saved atom");
      saved->SetPosition(-0.09, 0.00001, 12.9800007);
      saved->SetPartialCharge(0.000030006);
      saved->SetStereochemistry(AtomStereo::S);
      saved->SetAromaticity(true);
      bond->SetTag(4);
      angle->SetTag(5);
      dihedral->SetTag(6);
      fixture.atm.AddAngle(angle);
      fixture.atm.AddDihedral(dihedral);
      fixture.atm.AddBond(bond);
      saved->SetType(fixture.fftype);

      std::ostringstream os;
      {
        Out oar(os);
        check_nothrow(oar(saved, bond, angle, dihedral, fixture.mol,
    fixture.fftype));
      }

      Atom loaded;
      Bond bond_loaded;
      Angle angle_loaded;
      Dihedral dihedral_loaded;
      Molecule mol_loaded;
      FFAtom fftype_loaded;
      std::istringstream is(os.str());
      {
        In iar(is);
        check_nothrow(iar(loaded, bond_loaded, angle_loaded, dihedral_loaded,
                          mol_loaded, fftype_loaded));
      }

      fixture.atm.imp = loaded;
      check_eq(saved->GetElement(), loaded->GetElement());
      check_eq(saved->GetFormalCharge(), loaded->GetFormalCharge());
      check_eq(saved->GetTag(), loaded->GetTag());
      check_eq(saved->GetImplicitCount(), loaded->GetImplicitCount());
      check_eq(saved->GetName(), loaded->GetName());
      check_eq(saved->GetPosition(), loaded->GetPosition());
      check_eq(saved->GetPartialCharge(), loaded->GetPartialCharge());
      check_eq(saved->GetStereochemistry(), loaded->GetStereochemistry());
      check_eq(saved->GetAromaticity(), loaded->GetAromaticity());
      check_eq(bond->GetTag(), loaded->GetBondIters().first->lock()->GetTag());
      check_eq(angle->GetTag(),
    loaded->GetAngleIters().first->lock()->GetTag());
      check_eq(dihedral->GetTag(),
    loaded->GetDihedralIters().first->lock()->GetTag()); check_eq(bond_loaded,
    fixture.atm.get_bnds().front().lock()); check_eq(angle_loaded,
    fixture.atm.get_angs().front().lock()); check_eq(dihedral_loaded,
    fixture.atm.get_dhds().front().lock()); check_eq(fftype_loaded,
    loaded->GetType());
    }
    DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixatom_serial, ixserial<IXAtom>);
   */

  /*  test_case_fixture(test::AtomTestFixture, "IXAtom construction") {
      check_nothrow(test::TestAtom tatm(mol));
      test::TestAtom tatm(mol);

      check_eq(mol, tatm.get_mol().lock());
      check_eq(Element(), tatm.get_elem().lock());
      check_eq(0, tatm.get_fc());
      check_eq(0, tatm.get_tag());
      check_eq(0, tatm.get_implicitH());
      check_eq("", tatm.get_name());
      check_eq(Vec3(0.0,0.0,0.0), tatm.get_pos());
      check_eq(0.0, tatm.get_partial());
      check_eq(AtomStereo::UNDEFINED, tatm.get_stereo());
      check_eq(false, tatm.get_aromatic());
      check_eq(0, tatm.get_bnds().size());
      check_eq(0, tatm.get_angs().size());
      check_eq(0, tatm.get_dhds().size());
      check_eq(FFAtom(), tatm.get_type());

      // Check unique IDs correctly update
      test::TestAtom atm1(mol);
      test::TestAtom atm2(mol);
      check_ne(atm1.GetUniqueID(), atm2.GetUniqueID());
      check_eq(atm1.GetUniqueID() + 1, atm2.GetUniqueID());
    }
   */

  /*  test_case_fixture(test::AtomTestFixture, "IXAtom getting and setting") {
      // Check nothrow
      check_nothrow(atm.AddImplicitHydrogen());
      check_nothrow(atm.RemoveImplicitHydrogen());
      check_nothrow(atm.SetElement(GetPeriodicTable()->GetElement("Zr")));
      check_nothrow(atm.SetElement("Pb"));
      check_nothrow(atm.SetElement(32));
      check_nothrow(atm.SetFormalCharge(-2));
      check_nothrow(atm.SetPartialCharge(-0.09));
      check_nothrow(atm.SetImplicitCount(4));
      check_nothrow(atm.SetTag(2));
      check_nothrow(atm.SetName("TestAtom"));
      check_nothrow(atm.SetX(1.1));
      check_nothrow(atm.SetY(-2.2));
      check_nothrow(atm.SetZ(3.4));
      check_nothrow(atm.SetPosition(4.5, -7.9002, 0.00004));
      check_nothrow(atm.SetStereochemistry(AtomStereo::ACHIRAL));
      check_nothrow(atm.SetAromaticity(true));
      check_nothrow(atm.SetType(fftype));

      // Check correctness of gets
      check_eq(GetPeriodicTable()->GetElement(32), atm.GetElement());
      check_eq(GetPeriodicTable()->GetElement(32), atm.get_elem().lock());
      check_eq(-2, atm.GetFormalCharge());
      check_eq(-2, atm.get_fc());
      check_eq(approximately(-0.09), atm.GetPartialCharge());
      check_eq(approximately(-0.09), atm.get_partial());
      check_eq(2, atm.GetTag());
      check_eq(2, atm.get_tag());
      check_eq(4, atm.GetImplicitCount());
      check_eq(4, atm.get_implicitH());
      check_eq(mol, atm.GetMolecule());
      check_eq(mol, atm.get_mol().lock());
      check_eq("TestAtom", atm.GetName());
      check_eq("TestAtom", atm.get_name());
      check_eq(approximately(4.5), atm.GetX());
      check_eq(approximately(-7.9002), atm.GetY());
      check_eq(approximately(0.00004), atm.GetZ());
      check_eq(Vec3(4.5, -7.9002, 0.00004), atm.GetPosition());
      check_eq(Vec3(4.5, -7.9002, 0.00004), atm.get_pos());
      check_eq(AtomStereo::ACHIRAL, atm.GetStereochemistry());
      check_eq(AtomStereo::ACHIRAL, atm.get_stereo());
      check_eq(true, atm.GetAromaticity());
      check_eq(true, atm.get_aromatic());
      check_eq(atm.GetTag(), atm.GetIndex());
      check_eq(fftype, atm.GetType());

      // Check no owning of molecule
      mol.reset();
      check(atm.get_mol().expired());
      check_eq(Molecule(), atm.GetMolecule());

      // Check get index returns tag when molecule dead
      check_eq(atm.GetTag(), atm.GetIndex());

      // Check correctness of overloaded sets
      Element e1 = GetPeriodicTable()->GetElement(12);
      Element e2 = GetPeriodicTable()->GetElement(13);
      Element e3 = GetPeriodicTable()->GetElement(14);
      atm.SetElement(e1);
      check_eq(e1, atm.get_elem().lock());
      atm.SetElement(e2->GetName());
      check_eq(e2, atm.get_elem().lock());
      atm.SetElement(e3->GetAtomicNumber());
      check_eq(e3, atm.get_elem().lock());
    }
   */

  /*  test_case_fixture(test::AtomTestFixture, "IXAtom adding and removing") {
      //  Remember, normal use is only through IXMolecule, so some assumptions
      //  are made. These are that will never attempt to add the same item
    twice,
      //  and will never attempt to remove an item which has not been added.
      std::random_device rd;
      std::mt19937 generator(rd());
      std::uniform_int_distribution<size_t> distribution(7,17);
      size_t num = distribution(generator);

      subcase("implicit hydrogen") {
        size_t count = 0;
        while (atm.GetImplicitCount() < num)
          check_eq(++count, atm.AddImplicitHydrogen());
        while (atm.GetImplicitCount())
          check_eq(--count, atm.RemoveImplicitHydrogen());
        // Check removing implicit H when at 0 doesn't wrap around
        check_eq(0, atm.RemoveImplicitHydrogen());
      }

      subcase("bonds") {
        std::vector<Bond> bonds; bonds.reserve(num);
        // adding bonds
        for (size_t i = 1; i <= num; ++i) {
          bonds.push_back(test::CreateGenericTestBond().imp);
          check_nothrow(atm.AddBond(bonds.back()));
          check_eq(i, atm.NumBonds());
        }

        // checking state of added bonds
        check_eq(num, std::distance(atm.GetBondIters().first,
                                    atm.GetBondIters().second));
        check_eq(std::make_pair(atm.get_bnds().begin(), atm.get_bnds().end()),
                 atm.GetBondIters());
        std::vector<Bond> saved_bonds; saved_bonds.reserve(num);
        for (_Bond bnd : atm.get_bnds()) saved_bonds.emplace_back(bnd.lock());
        check_eq(bonds.size(), saved_bonds.size());
        check_eq(bonds, saved_bonds);

        // removing bonds (randomise removal order)
        std::shuffle(bonds.begin(), bonds.end(), generator);
        for (Bond bnd : bonds) {
          check_nothrow(atm.RemoveBond(bnd));
          check_eq(--num, atm.NumBonds());
          check_eq(num, std::distance(atm.get_bnds().begin(),
                                      atm.get_bnds().end()));
          check_eq(std::make_pair(atm.get_bnds().begin(), atm.get_bnds().end()),
                   atm.GetBondIters());
        }
      }

      subcase("angles") {
        std::vector<Angle> angles; angles.reserve(num);
        for (size_t i = 1; i <= num; ++i) {
          angles.push_back(test::CreateGenericTestAngle().imp);
          check_nothrow(atm.AddAngle(angles.back()));
          check_eq(i, atm.NumAngles());
        }

        // checking state of added angles
        check_eq(num, std::distance(atm.GetAngleIters().first,
                                    atm.GetAngleIters().second));
        check_eq(std::make_pair(atm.get_angs().begin(), atm.get_angs().end()),
                 atm.GetAngleIters());
        std::vector<Angle> saved_angles; saved_angles.reserve(num);
        for (_Angle ang : atm.get_angs()) saved_angles.emplace_back(ang.lock());
        check_eq(angles.size(), saved_angles.size());
        check_eq(angles, saved_angles);

        // removing angles (randomise removal order)
        std::shuffle(angles.begin(), angles.end(), generator);
        for (Angle ang : angles) {
          check_nothrow(atm.RemoveAngle(ang));
          check_eq(--num, atm.NumAngles());
          check_eq(num, std::distance(atm.get_angs().begin(),
                                      atm.get_angs().end()));
          check_eq(std::make_pair(atm.get_angs().begin(), atm.get_angs().end()),
                   atm.GetAngleIters());
        }
      }

      subcase("dihedrals") {
        std::vector<Dihedral> dihedrals; dihedrals.reserve(num);
        for (size_t i = 1; i <= num; ++i) {
          dihedrals.push_back(test::CreateGenericTestDihedral().imp);
          check_nothrow(atm.AddDihedral(dihedrals.back()));
          check_eq(i, atm.NumDihedrals());
        }

        // checking state of added dihedrals
        check_eq(num, std::distance(atm.GetDihedralIters().first,
                                    atm.GetDihedralIters().second));
        check_eq(std::make_pair(atm.get_dhds().begin(), atm.get_dhds().end()),
                 atm.GetDihedralIters());
        std::vector<Dihedral> saved_dihedrals; saved_dihedrals.reserve(num);
        for (_Dihedral dhd : atm.get_dhds())
    saved_dihedrals.emplace_back(dhd.lock()); check_eq(dihedrals.size(),
    saved_dihedrals.size()); check_eq(dihedrals, saved_dihedrals);

        // removing dihedrals (randomise removal order)
        std::shuffle(dihedrals.begin(), dihedrals.end(), generator);
        for (Dihedral dhd : dihedrals) {
          check_nothrow(atm.RemoveDihedral(dhd));
          check_eq(--num, atm.NumDihedrals());
          check_eq(num, std::distance(atm.get_dhds().begin(),
                                      atm.get_dhds().end()));
          check_eq(std::make_pair(atm.get_dhds().begin(), atm.get_dhds().end()),
                   atm.GetDihedralIters());
        }
      }

    }
   */

  /*  test_case_fixture(test::AtomTestFixture, "IXAtom printing methods") {
      std::stringstream ss;
      check_nothrow(ss << atm.imp);
      check_eq("Atom(0)", ss.str());
      check_eq("Atom(0, XX)", atm.ToString());
      ss.str("");

      atm.SetTag(12); atm.SetElement("Zr");
      check_nothrow(ss << atm.imp);
      check_eq("Atom(12)", ss.str());
      check_eq("Atom(12, Zr)", atm.ToString());
      ss.str("");

      // Check blank atom doesn't do anything
      check_nothrow(ss << Atom());
      check_eq("", ss.str());
    }
    */

  test_suite_close();
} // namespace indigox
