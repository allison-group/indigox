#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/quad.hpp>
#include <indigox/utils/serialise.hpp>
#include <indigox/utils/triple.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/atom_test.hpp>
#include <indigox/test/angle_test.hpp>
#include <indigox/test/bond_test.hpp>
#include <indigox/test/dihedral_test.hpp>
#include <indigox/test/molecule_test.hpp>

namespace indigox {
  test_suite_open("Molecule");
  
// ============================================================================
// == SERIALISATION ===========================================================
// ============================================================================
  
  template <typename Archive>
  void Molecule::save(Archive& archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("atoms", _atms),
            INDIGOX_SERIAL_NVP("bonds", _bnds),
            INDIGOX_SERIAL_NVP("angles", _angs),
            INDIGOX_SERIAL_NVP("dihedrals", _dhds),
            INDIGOX_SERIAL_NVP("name", _name),
            INDIGOX_SERIAL_NVP("molecular_charge", _q),
            INDIGOX_SERIAL_NVP("molecular_graph", _g));
  }
  
  template <typename Archive>
  void Molecule::load(Archive& archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("atoms", _atms),
            INDIGOX_SERIAL_NVP("bonds", _bnds),
            INDIGOX_SERIAL_NVP("angles", _angs),
            INDIGOX_SERIAL_NVP("dihedrals", _dhds),
            INDIGOX_SERIAL_NVP("name", _name),
            INDIGOX_SERIAL_NVP("molecular_charge", _q),
            INDIGOX_SERIAL_NVP("molecular_graph", _g));
  }
  INDIGOX_SERIALISE_SPLIT(Molecule);
  
/*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMolecule serialisation", T, ixmol_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::MoleculeTestFixture fixture;
    fixture.BuildBenzene();
    fixture.benzene.PerceiveAngles();
    fixture.benzene.PerceiveDihedrals();
    Molecule saved = fixture.benzene.imp;
    saved->SetName("Benzene");
    saved->SetMolecularCharge(12);

    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved));
    }

    Molecule loaded;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded));
    }
    fixture.blankmol.imp = loaded;

    // loaded system should have emergent all set, and no cached values
    // These checks come first as NumAngles and NumDihedrals call their perceive
    // methods, resetting the values.
    check_eq("", fixture.blankmol.get_formula_cache());
    test::MoleculeTestFixture::EmergeSet eset; eset.set();
    check_eq(eset, fixture.blankmol.get_emerge());

    // Check contains are right
    check_eq(saved->GetName(), loaded->GetName());
    check_eq(saved->GetMolecularCharge(), loaded->GetMolecularCharge());
    check_eq(saved->NumAtoms(), loaded->NumAtoms());
    check_eq(saved->NumBonds(), loaded->NumBonds());
    check_eq(saved->NumAngles(), loaded->NumAngles());
    check_eq(saved->NumDihedrals(), loaded->NumDihedrals());
    check(loaded->GetGraph());
    check_eq(saved->GetGraph()->NumVertices(), loaded->GetGraph()->NumVertices());
    check_eq(saved->GetGraph()->NumEdges(), loaded->GetGraph()->NumEdges());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixmol_serial, ixserial<IXMolecule>);
 */
  
// ============================================================================
// == CONSTRUCTION AND INITALISATION ==========================================
// ============================================================================
  
  Molecule::Molecule() :
  utils::IXCountableObject<Molecule>(), utils::ModifiableObject(), _name(""),
  _q(0), _g(std::make_shared<graph::MolecularGraph>()), _formula_cache(0, ""),
  _angle_perceive(0), _dihedral_perceive(0) { }
  
  void Molecule::Init() {
    _g->_mol = weak_from_this();
  }
  
  Molecule::~Molecule() {
    _g->Clear();
    for (sAtom atm : _atms) atm->Clear();
    _atms.clear();
    for (sBond bnd : _bnds) bnd->Clear();
    _bnds.clear();
    for (sAngle ang : _angs) ang->Clear();
    _angs.clear();
    for (sDihedral dhd : _dhds) dhd->Clear();
    _dhds.clear();
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule construction") {
    check_nothrow(test::TestMolecule tmol);

    check_eq("", blankmol.get_name());
    check_eq(0, blankmol.get_q());
    check_eq(0, blankmol.get_atms().size());
    check_eq(0, blankmol.get_bnds().size());
    check_eq(0, blankmol.get_angs().size());
    check_eq(0, blankmol.get_dhds().size());
    check_eq(test::TestMolecule::EmergeSet().flip(), blankmol.get_emerge());
    check_eq(graph::MolecularGraph(), blankmol.get_g());
    check_eq("", blankmol.get_formula_cache());

    // Check init does what it should
    check_nothrow(blankmol.Init());
    check_ne(graph::MolecularGraph(), blankmol.get_g());

    // Check unique IDs correctly update
    test::TestMolecule mol1, mol2;
    check_ne(mol1.GetUniqueID(), mol2.GetUniqueID());
    check_eq(mol1.GetUniqueID() + 1, mol2.GetUniqueID());
  }
 */
  
// ============================================================================
// == FINDING ITEMS ===========================================================
// ============================================================================
  
  int64_t Molecule::_FindBond(Atom &a, Atom &b) const {
    auto pred = [&a, &b](sBond bnd) {
      auto s = bnd->GetAtoms();
      return ((&s.first == &a && &s.second == &b)
              || (&s.first == &b && &s.second == &a));
    };
    MolBondIter pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return pos == _bnds.end() ? -1 : std::distance(_bnds.begin(), pos);
  }
  
  int64_t Molecule::_FindAngle(Atom &a, Atom &b, Atom &c) const {
    auto pred = [&a, &b, &c] (sAngle ang) {
      auto s = ang->GetAtoms();
      return (&s.second == &b) && ((&s.first == &a  && &s.third == &c)
                                   || (&s.first == &c  && &s.third == &a));
    };
    MolAngleIter pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return pos == _angs.end() ? -1 : std::distance(_angs.begin(), pos);
  }
  
  int64_t Molecule::_FindDihedral(Atom &a, Atom &b, Atom &c, Atom &d) const {
    auto pred = [&a, &b, &c, &d] (sDihedral dhd) {
      auto s = dhd->GetAtoms();
      return ((&s.first == &a && &s.second == &b && &s.third == &c && &s.fourth == &d)
              || (&s.first == &d && &s.second == &c && &s.third == &b && &s.fourth == &a));
    };
    MolDihedralIter pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    return pos == _dhds.end() ? -1 : std::distance(_dhds.begin(), pos);
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule find items") {
    BuildBenzene();
    BuildRandomMolecule();
    subcase("finding bonds") {
      std::shuffle(randmol.get_bnds().begin(), randmol.get_bnds().end(),
                   generator);  // mix up the bonds order
      for (size_t i = 0; i < b_randmol.size(); ++i) {
        Atom a = a_randmol[e_randmol[i].first];
        Atom b = a_randmol[e_randmol[i].second];
        check_eq(b_randmol[i], randmol.FindBond(a, b));
      }
      // Finding a bond that doesn't exist should return an empty bond
      check_eq(Bond(), randmol.FindBond(a_benzene[0], a_benzene[1]));
      check_eq(Bond(), randmol.FindBond(a_benzene[2], Atom()));
      check_eq(Bond(), randmol.FindBond(Atom(), a_benzene[3]));
    }

    subcase("finding angles") {
      benzene.PerceiveAngles();
      std::set<Angle> angles(benzene.GetAngles().first,
                             benzene.GetAngles().second);
      std::set<Angle> got_angles;
      std::shuffle(benzene.get_angs().begin(), benzene.get_angs().end(),
                   generator);
      for (size_t i = 0; i < g_benzene.size(); ++i) {
        Atom a = a_benzene[g_benzene[i].first];
        Atom b = a_benzene[g_benzene[i].second];
        Atom c = a_benzene[g_benzene[i].third];
        Angle ang = benzene.FindAngle(a, b, c);
        check_ne(Angle(), ang);
        got_angles.insert(ang);
      }
      check_eq(angles, got_angles);
      check_eq(got_angles.size(), g_benzene.size());
      // Finding an angle that doesn't exist should return an empty angle
      check_eq(Angle(), randmol.FindAngle(a_randmol[0], a_randmol[3],
                                          a_randmol[2]));
      check_eq(Angle(), randmol.FindAngle(Atom(), a_randmol[4],
                                          a_randmol[5]));
      check_eq(Angle(), randmol.FindAngle(a_randmol[6], Atom(),
                                          a_randmol[7]));
    }

    subcase("finding dihedrals") {
      benzene.PerceiveDihedrals();
      std::set<Dihedral> dihedrals(benzene.GetDihedrals().first,
                                   benzene.GetDihedrals().second);
      std::set<Dihedral> got_dihedrals;
      std::shuffle(benzene.get_dhds().begin(), benzene.get_dhds().end(),
                   generator);
      for (size_t i = 0; i < d_benzene.size(); ++i) {
        Atom a = a_benzene[d_benzene[i].first];
        Atom b = a_benzene[d_benzene[i].second];
        Atom c = a_benzene[d_benzene[i].third];
        Atom d = a_benzene[d_benzene[i].fourth];
        Dihedral dhd = benzene.FindDihedral(a, b, c, d);
        check_ne(Dihedral(), dhd);
        got_dihedrals.insert(dhd);
      }
      check_eq(dihedrals, got_dihedrals);
      check_eq(got_dihedrals.size(), d_benzene.size());
      // Finding a dihedral that doesn't exist should return an empty dihedral
      check_eq(Dihedral(), randmol.FindDihedral(a_randmol[0],
                                                a_randmol[1],
                                                a_randmol[2],
                                                a_randmol[3]));
      check_eq(Dihedral(), randmol.FindDihedral(a_randmol[0], Atom(),
                                                a_randmol[2],
                                                a_randmol[3]));
      check_eq(Dihedral(), randmol.FindDihedral(a_randmol[0],
                                                a_randmol[1],
                                                a_randmol[2], Atom()));
    }
  }
  */
  
  bool Molecule::HasBond(Atom& a, Atom& b) const {
    if (!HasAtom(a) || !HasAtom(b)) return false;
    return _FindBond(a, b) != -1;
  }
  
  bool Molecule::HasAngle(Atom &a, Atom &b, Atom &c) {
    if (!HasAtom(a) || !HasAtom(b) || !HasAtom(c)) return false;
    PerceiveAngles();
    return _FindAngle(a, b, c) != -1;
  }
  
  bool Molecule::HasDihedral(Atom &a, Atom &b, Atom &c, Atom &d) {
    if (!HasAtom(a) || !HasAtom(b) || !HasAtom(c) || !HasAtom(d)) return false;
    PerceiveDihedrals();
    return _FindDihedral(a, b, c, d) != -1;
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule has items") {
    BuildBenzene();
    BuildRandomMolecule();

    subcase("has atoms") {
      for (Atom atm : a_randmol) {
        check_false(benzene.HasAtom(atm));
        check(randmol.HasAtom(atm));
      }
      // Null check
      check_false(benzene.HasAtom(Atom()));
    }

    subcase("has bonds") {
      for (Bond bnd : b_randmol) {
        check_false(benzene.HasBond(bnd));
        check(randmol.HasBond(bnd));
      }

      size_t i, j;
      for (std::pair<size_t,size_t> bnd : e_randmol) {
        std::tie(i, j) = bnd;
        check(randmol.HasBond(a_randmol[i], a_randmol[j]));
        check_false(benzene.HasBond(a_randmol[i], a_randmol[j]));
      }
      check_false(benzene.HasBond(a_benzene[6], a_benzene[7]));
      check_false(benzene.HasBond(Bond()));
    }

    subcase("has angles") {
      // Perceive angles automatically when check for angle by atoms
      check_eq(0, benzene.get_angs().size());
      check_false(benzene.HasAngle(a_benzene[6],
                                   a_benzene[7],
                                   a_benzene[8]));
      check_eq(g_benzene.size(), benzene.get_angs().size());
      check_false(randmol.HasAngle(Angle()));
      randmol.PerceiveAngles();

      for (Angle ang : randmol.get_angs()) {
        check(randmol.HasAngle(ang));
        check_false(benzene.HasAngle(ang));
      }

      for (auto ang : g_benzene) {
        check(benzene.HasAngle(a_benzene[ang.first],
                               a_benzene[ang.second],
                               a_benzene[ang.third]));
        check_false(randmol.HasAngle(a_benzene[ang.first],
                                     a_benzene[ang.second],
                                     a_benzene[ang.third]));
      }
      check_false(benzene.HasAngle(a_benzene[6],
                                   a_benzene[7],
                                   a_benzene[8]));
    }

    subcase("has dihedrals") {
      // Perceive dihedrals automatically when check for dihedral by atoms
      check_eq(0, benzene.get_dhds().size());
      check_false(benzene.HasDihedral(a_benzene[6], a_benzene[7],
                                      a_benzene[8], a_benzene[9]));
      check_eq(d_benzene.size(), benzene.get_dhds().size());
      check_false(randmol.HasDihedral(Dihedral()));
      randmol.PerceiveDihedrals();
      for (Dihedral dhd : randmol.get_dhds()) {
        check(randmol.HasDihedral(dhd));
        check_false(benzene.HasDihedral(dhd));
      }

      for (auto dhd : d_benzene) {
        check(benzene.HasDihedral(a_benzene[dhd.first],
                                  a_benzene[dhd.second],
                                  a_benzene[dhd.third],
                                  a_benzene[dhd.fourth]));
        check_false(randmol.HasDihedral(a_benzene[dhd.first],
                                        a_benzene[dhd.second],
                                        a_benzene[dhd.third],
                                        a_benzene[dhd.fourth]));
      }
    }
  }
*/
  
// ============================================================================
// == GETTING AND SETTING =====================================================
// ============================================================================
  
  Atom& Molecule::GetAtomTag(uint32_t tag) const {
    auto pred = [tag](sAtom a) { return a->GetTag() == tag; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    if (pos == _atms.end()) throw std::out_of_range("No atom with tag found");
    return *(*pos);
  }
  
  Atom& Molecule::GetAtomID(uint32_t id) const {
    auto pred = [id](sAtom a) { return a->GetUniqueID() == id; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    if (pos == _atms.end()) throw std::out_of_range("No atom with ID found");
    return *(*pos);
  }
  
  Bond& Molecule::GetBond(Atom& a, Atom& b) const {
    int64_t pos = _FindBond(a, b);
    if (pos == -1) throw std::out_of_range("No bond between atoms");
    return *_bnds[pos];
  }
  
  Bond& Molecule::GetBondTag(uint32_t tag) const {
    auto pred = [tag](sBond b) { return b->GetTag() == tag; };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    if (pos == _bnds.end()) throw std::out_of_range("No bond with tag found");
    return *(*pos);
  }
  
  Bond& Molecule::GetBondID(uint32_t id) const {
    auto pred = [id](sBond b) { return b->GetUniqueID() == id; };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    if (pos == _bnds.end()) throw std::out_of_range("No bond with ID found");
    return *(*pos);
  }
  
  Angle& Molecule::GetAngle(Atom& a, Atom&b, Atom& c) {
    PerceiveAngles();
    int64_t pos = _FindAngle(a, b, c);
    if (pos == -1) throw std::out_of_range("No angle between atoms");
    return *_angs[pos];
  }
  
  Angle& Molecule::GetAngleTag(uint32_t tag) const {
    auto pred = [tag](sAngle ang) { return ang->GetTag() == tag; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    if (pos == _angs.end()) throw std::out_of_range("No angle with tag found");
    return *(*pos);
  }
  
  Dihedral& Molecule::GetDihedral(Atom &a, Atom &b, Atom &c, Atom &d) {
    PerceiveDihedrals();
    int64_t pos = _FindDihedral(a, b, c, d);
    if (pos == -1) throw std::out_of_range("No dihedral between atoms");
    return *_dhds[pos];
  }
  
  Dihedral& Molecule::GetDihedralTag(uint32_t tag) const {
    auto pred = [tag](sDihedral dhd) { return dhd->GetTag() == tag; };
    auto pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    if (pos == _dhds.end()) throw std::out_of_range("No dihedral with tag found");
    return *(*pos);
  }
  
  Angle& Molecule::GetAngleID(uint32_t id) const {
    auto pred = [id](sAngle ang) { return ang->GetUniqueID() == id; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    if (pos == _angs.end()) throw std::out_of_range("No angle with ID found");
    return *(*pos);
  }
  
  Dihedral& Molecule::GetDihedralID(uint32_t id) const {
    auto pred = [id](sDihedral dhd) { return dhd->GetUniqueID() == id; };
    auto pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    if (pos == _dhds.end()) throw std::out_of_range("No dihedral with ID found");
    return *(*pos);
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule getting and setting") {
    BuildBenzene();
    benzene.PerceiveAngles();
    benzene.PerceiveDihedrals();

    benzene.get_emerge().reset();
    test::TestMolecule::EmergeSet empty_emergent_state;
    empty_emergent_state.reset();
    // Setting name should not change emergent state
    check_nothrow(benzene.SetName("Benzene"));
    check_eq("Benzene", benzene.GetName());
    check_eq(empty_emergent_state, benzene.get_emerge());

    // Setting molecular charge should change emergent state
    benzene.get_emerge().reset();
    check_nothrow(benzene.SetMolecularCharge(-5));
    check_eq(-5, benzene.GetMolecularCharge());
    check_ne(empty_emergent_state, benzene.get_emerge());
    // But only when the value actually changes
    benzene.get_emerge().reset();
    benzene.SetMolecularCharge(-5);
    check_eq(empty_emergent_state, benzene.get_emerge());

    // Getting graph should not be empty
    check_ne(graph::MolecularGraph(), benzene.GetGraph());
    check_eq(benzene.NumAtoms(), benzene.GetGraph()->NumVertices());
    check_eq(benzene.NumBonds(), benzene.GetGraph()->NumEdges());

    // -- ATOMS --
    // Getting atom by position
    for (size_t i = 0; i < a_benzene.size(); ++i)
      check_eq(a_benzene[i], benzene.GetAtom(i));
    // Over position should return null atom
    check_eq(Atom(), benzene.GetAtom(a_benzene.size()));

    // Getting atom by ID/tag
    for (Atom atm : a_benzene) {
      check_eq(atm, benzene.GetAtomID(atm->GetUniqueID()));
      check_eq(atm, benzene.GetAtomTag(atm->GetTag()));
    }
    // Bad ID/tag should return null
    check_eq(Atom(), benzene.GetAtomID(std::numeric_limits<uint32_t>::max()));
    check_eq(Atom(), benzene.GetAtomTag(std::numeric_limits<uint32_t>::max()));
    // Should return first atom with given tag
    a_benzene.front()->SetTag(a_benzene.back()->GetTag());
    check_eq(a_benzene.front(), benzene.GetAtomTag(a_benzene.back()->GetTag()));

    // Getting atoms by iterator
    auto atm_iters = benzene.GetAtoms();
    std::vector<Atom> got_atoms(atm_iters.first, atm_iters.second);
    check_eq(a_benzene, got_atoms);

    // -- BONDS --
    // Getting bond by position/ atom pairs
    for (size_t i = 0; i < b_benzene.size(); ++i) {
      check_eq(b_benzene[i], benzene.GetBond(i));
      auto pos = e_benzene[i];
      check_eq(b_benzene[i], benzene.GetBond(a_benzene[pos.first],
                                             a_benzene[pos.second]));
    }
    // Over position should return null bond
    check_eq(Bond(), benzene.GetBond(a_benzene.size()));

    // Getting bond by ID/tag
    for (Bond bnd : b_benzene) {
      check_eq(bnd, benzene.GetBondID(bnd->GetUniqueID()));
      check_eq(bnd, benzene.GetBondTag(bnd->GetTag()));
    }
    // Bad ID/Tag should return null
    check_eq(Bond(), benzene.GetBondID(std::numeric_limits<uint32_t>::max()));
    check_eq(Bond(), benzene.GetBondTag(std::numeric_limits<uint32_t>::max()));
    // Should return first bond with given tag
    b_benzene.front()->SetTag(b_benzene.back()->GetTag());
    check_eq(b_benzene.front(), benzene.GetBondTag(b_benzene.back()->GetTag()));

    // Getting bonds by iterator
    auto bnd_iters = benzene.GetBonds();
    std::vector<Bond> got_bonds(bnd_iters.first, bnd_iters.second);
    check_eq(b_benzene, got_bonds);

    // -- ANGLES --
    // Getting angle by position/atom triple
    std::vector<Angle> p_angles(benzene.get_angs().begin(),
                                benzene.get_angs().end());
    for (size_t i = 0; i < g_benzene.size(); ++i) {
      check_eq(p_angles[i], benzene.GetAngle(i));
      auto pos = g_benzene[i];
      check_eq(p_angles[i], benzene.GetAngle(a_benzene[pos.third],
                                             a_benzene[pos.second],
                                             a_benzene[pos.first]));
    }
    // Over position should return null angle
    check_eq(Angle(), benzene.GetAngle(g_benzene.size()));

    // Getting angle by ID/tag
    size_t new_tag = p_angles.size() + 3;
    for (Angle ang : p_angles) {
      ang->SetTag(++new_tag);
      check_eq(ang, benzene.GetAngleID(ang->GetUniqueID()));
      check_eq(ang, benzene.GetAngleTag(ang->GetTag()));
    }
    // Bad ID/tag should return null
    check_eq(Angle(), benzene.GetAngleID(std::numeric_limits<uint32_t>::max()));
    check_eq(Angle(), benzene.GetAngleTag(std::numeric_limits<uint32_t>::max()));
    // Should return first angle with given tag
    p_angles.front()->SetTag(p_angles.back()->GetTag());
    check_eq(p_angles.front(), benzene.GetAngleTag(p_angles.back()->GetTag()));

    // Getting angles by iterator
    auto ang_iters = benzene.GetAngles();
    std::vector<Angle> got_angles(ang_iters.first, ang_iters.second);
    check_eq(p_angles, got_angles);

    // -- DIHEDRALS --
    // Getting dihedral by position/atom quad
    std::vector<Dihedral> p_dihedrals(benzene.get_dhds().begin(),
                                      benzene.get_dhds().end());
    for (size_t i = 0; i < d_benzene.size(); ++i) {
      check_eq(p_dihedrals[i], benzene.GetDihedral(i));
      auto pos = d_benzene[i];
      check_eq(p_dihedrals[i], benzene.GetDihedral(a_benzene[pos.fourth],
                                                   a_benzene[pos.third],
                                                   a_benzene[pos.second],
                                                   a_benzene[pos.first]));
    }
    // Over position should return null dihedral
    check_eq(Dihedral(), benzene.GetDihedral(d_benzene.size()));

    // Getting dihedral by ID/tag
    for (Dihedral dhd : p_dihedrals) {
      dhd->SetTag(++new_tag);
      check_eq(dhd, benzene.GetDihedralID(dhd->GetUniqueID()));
      check_eq(dhd, benzene.GetDihedralTag(dhd->GetTag()));
    }
    // Bad ID/tag should return null
    check_eq(Dihedral(), benzene.GetDihedralID(std::numeric_limits<uint32_t>::max()));
    check_eq(Dihedral(), benzene.GetDihedralTag(std::numeric_limits<uint32_t>::max()));
    // SHould return first dihedral with given tag
    p_dihedrals.front()->SetTag(new_tag);
    check_eq(p_dihedrals.front(), benzene.GetDihedralTag(new_tag));

    // Getting dihedrals by iterator
    auto dhd_iters = benzene.GetDihedrals();
    std::vector<Dihedral> got_dihedrals(dhd_iters.first, dhd_iters.second);
    check_eq(p_dihedrals, got_dihedrals);

    // Numbers of things should be right
    check_eq(12, benzene.NumAtoms());
    check_eq(12, benzene.NumBonds());
    check_eq(18, benzene.NumAngles());
    check_eq(24, benzene.NumDihedrals());
  }
*/
  
  std::string Molecule::GetFormula() {
    ModifiableObject::State state = GetCurrentState();
    if (state != _formula_cache.first || _formula_cache.first == 0) {
      std::map<std::string, size_t> e_count;
      for (sAtom atm : _atms) e_count[atm->GetElement().GetSymbol()]++;
      std::stringstream ss;
      if (e_count["C"]) ss << "C";
      if (e_count["C"] > 1) ss << e_count["C"];
      if (e_count["H"]) ss << "H";
      if (e_count["H"] > 1) ss << e_count["H"];
      for (auto& e : e_count) {
        if (e.first != "C" && e.first != "H") {
          ss << e.first;
          if (e.second > 1) ss << e.second;
        }
      }
      _formula_cache = std::make_pair(state, ss.str());
    }
    return _formula_cache.second;
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule getting formula") {
    BuildBenzene();
    size_t e_pos = static_cast<size_t>(EmergeProp::MOLECULAR_FORMULA);
    
    // Check correct calculation and setting of emergent property state
    check(benzene.get_emerge()[e_pos]);
    check_eq("C6H2BrClFI", benzene.GetFormula());
    check_false(benzene.get_emerge()[e_pos]);
    
    // Check formula is cached
    check_eq("C6H2BrClFI", benzene.get_formula_cache());
    
    // Check recalculates only when emergent property set
    a_benzene.front()->SetElement("F");
    benzene.get_emerge().reset(e_pos);
    check_eq("C6H2BrClFI", benzene.GetFormula()); // shouldn't recalculate
    benzene.get_emerge().set(e_pos);
    check_eq("C5H2BrClF2I", benzene.GetFormula()); // should recalculate
  }
  
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule get part index") {
    BuildRandomMolecule();
    std::shuffle(randmol.get_atms().begin(), randmol.get_atms().end(), generator);
    for (size_t i = 0; i < a_randmol.size(); ++i) {
      size_t pos = randmol.get_atms()[i]->GetIndex();
      check_ne(pos, randmol.get_atms()[i]->GetTag());
      check_eq(i, pos);
    }
    
    std::shuffle(randmol.get_bnds().begin(), randmol.get_bnds().end(), generator);
    for (size_t i = 0; i < b_randmol.size(); ++i) {
      size_t pos = randmol.get_bnds()[i]->GetIndex();
      check_ne(pos, randmol.get_bnds()[i]->GetTag());
      check_eq(i, pos);
    }
    
    randmol.PerceiveAngles();
    std::shuffle(randmol.get_angs().begin(), randmol.get_angs().end(), generator);
    for (size_t i = 0; i < randmol.NumAngles(); ++i)
      check_eq(i, randmol.get_angs()[i]->GetIndex());
    
    randmol.PerceiveDihedrals();
    std::shuffle(randmol.get_dhds().begin(), randmol.get_dhds().end(), generator);
    for (size_t i = 0; i < randmol.NumDihedrals(); ++i)
      check_eq(i, randmol.get_dhds()[i]->GetIndex());
  }
 */
 
// ============================================================================
// == PERCEPTION ==============================================================
// ============================================================================
  
  size_t Molecule::PerceiveAngles() {
    ModifiableObject::State state = GetCurrentState();
    if (state && state == _angle_perceive) return 0;
    
    // Expected number of angles
    auto sum = [&](size_t current, sAtom v) -> size_t {
      size_t degree = v->NumBonds();
      if (degree < 2) return current;
      return current + degree * (degree - 1) / 2;
    };
    size_t count = std::accumulate(_atms.begin(), _atms.end(), 0, sum);
    _angs.reserve(count);
    count = 0;
    
    std::vector<Atom*> nbrs;
    nbrs.reserve(10);
    // Adding new angles
    for (sAtom at : _atms) {
      if (at->NumBonds() < 2) continue;
      for (wBond bn : at->GetBonds()) {
        auto atms = bn.lock()->GetAtoms();
        if (&atms.first == at.get()) nbrs.emplace_back(&atms.second);
        else nbrs.emplace_back(&atms.first);
      }
      for (size_t i = 0; i < nbrs.size() - 1; ++i) {
        for (size_t j = i + 1; j < nbrs.size(); ++j) {
          if (_FindAngle(*nbrs[i], *at, *nbrs[j]) != -1) continue;
          NewAngle(*nbrs[i], *at, *nbrs[j]);
          ++count;
        }
      }
      nbrs.clear();
    }
    _angle_perceive = state;
    return count;
  }
  
  size_t Molecule::PerceiveDihedrals() {
    ModifiableObject::State state = GetCurrentState();
    if (state && state == _dihedral_perceive) return 0;
    
    // Expected number of dihedrals
    auto sum = [&](size_t current, sBond b) -> size_t {
      size_t b_degree = b->GetSourceAtom().NumBonds();
      if (b_degree < 2) return current;
      size_t c_degree = b->GetTargetAtom().NumBonds();
      if (c_degree < 2) return current;
      return current + (b_degree - 1) * (c_degree - 1);
    };
    size_t count = std::accumulate(_bnds.begin(), _bnds.end(), 0, sum);
    _dhds.reserve(count);
    count = 0;
    
    // Adding new dihedrals
    std::vector<Atom*> B_nbrs, C_nbrs;
    B_nbrs.reserve(10); C_nbrs.reserve(10);
    for (sBond bn : _bnds) {
      Atom& B = bn->GetAtoms().first;
      Atom& C = bn->GetAtoms().second;
      if (B.NumBonds() < 2 || B.NumBonds() < 2) continue;
      
      for (wBond bn : B.GetBonds()) {
        auto atms = bn.lock()->GetAtoms();
        if (&atms.first == &B) B_nbrs.emplace_back(&atms.second);
        else B_nbrs.emplace_back(&atms.first);
      }
      for (wBond bn : C.GetBonds()) {
        auto atms = bn.lock()->GetAtoms();
        if (&atms.first == &C) C_nbrs.emplace_back(&atms.second);
        else C_nbrs.emplace_back(&atms.first);
      }
      
      for (size_t i = 0; i < B_nbrs.size(); ++i) {
        if (B_nbrs[i] == &C) continue;
        for (size_t j = 0; j < C_nbrs.size(); ++j) {
          if (C_nbrs[j] == &B) continue;
          if (_FindDihedral(*B_nbrs[i], B, C, *C_nbrs[j]) != -1) continue;
          NewDihedral(*B_nbrs[i], B, C, *C_nbrs[j]);
          ++count;
        }
      }
      B_nbrs.clear(); C_nbrs.clear();
    }
    _dihedral_perceive = state;
    return count;
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule perception") {
    BuildBenzene();
    subcase("Perception of angles") {
      // Should add x angles
      check_eq(g_benzene.size(), benzene.PerceiveAngles());
      for (auto abc : g_benzene) {
        Atom a = a_benzene[abc.first];
        Atom b = a_benzene[abc.second];
        Atom c = a_benzene[abc.third];
        // All expected angles should exist
        check(benzene.HasAngle(a, b, c));
      }
      // Repercieving shouldnt add any more
      check_eq(0, benzene.PerceiveAngles());
      
      // Adding an atom should mean one more angle
      Atom added = benzene.NewAtom();
      benzene.NewBond(added, a_benzene[11]);
      subcase("direct reperception") {
        check_eq(1, benzene.PerceiveAngles());
        check(benzene.HasAngle(added, a_benzene[11], a_benzene[5]));
      }
      subcase("reperception through check") {
        check(benzene.HasAngle(added, a_benzene[11], a_benzene[5]));
      }
      subcase("reperception through num angles") {
        check_eq(g_benzene.size() + 1, benzene.NumAngles());
      }
    }
    
    subcase("Perception of dihedrals") {
      check_eq(d_benzene.size(), benzene.PerceiveDihedrals());
      for (auto abcd : d_benzene) {
        Atom a = a_benzene[abcd.first];
        Atom b = a_benzene[abcd.second];
        Atom c = a_benzene[abcd.third];
        Atom d = a_benzene[abcd.fourth];
        // All expected dihedrals should exist
        check(benzene.HasDihedral(a, b, c, d));
      }
      // Reperceiving shouldn't add anymore dihedrals
      check_eq(0, benzene.PerceiveDihedrals());
      
      // Adding an atom should mean two more dihedrals
      Atom added = benzene.NewAtom();
      benzene.NewBond(added, a_benzene[11]);
      subcase("direct reperception") {
        check_eq(2, benzene.PerceiveDihedrals());
        check(benzene.HasDihedral(added, a_benzene[11],
                                  a_benzene[5], a_benzene[0]));
        check(benzene.HasDihedral(added, a_benzene[11],
                                  a_benzene[5], a_benzene[4]));
      }
      subcase("reperception through check") {
        check(benzene.HasDihedral(added, a_benzene[11],
                                  a_benzene[5], a_benzene[0]));
        check(benzene.HasDihedral(added, a_benzene[11],
                                  a_benzene[5], a_benzene[4]));
      }
      subcase("reperception through num angles") {
        check_eq(d_benzene.size() + 2, benzene.NumDihedrals());
      }
    }
  }
 */
 
// ============================================================================
// == STRUCTURE MODIFICATION ==================================================
// ============================================================================
  
  Atom& Molecule::NewAtom() {
    sAtom atom = std::make_shared<Atom>(*this);
    _atms.emplace_back(atom);
    _g->AddVertex(*atom);
    ModificationMade();
    return *atom;
  }
  
  Atom& Molecule::NewAtom(Element& element) {
    Atom& atom = NewAtom();
    atom.SetElement(element);
    return atom;
  }
  
  Atom& Molecule::NewAtom(std::string name) {
    Atom& atom = NewAtom();
    atom.SetName(name);
    return atom;
  }
  
  Atom& Molecule::NewAtom(std::string name, Element& element) {
    Atom& atom = NewAtom();
    atom.SetName(name);
    atom.SetElement(element);
    return atom;
  }
  
  bool Molecule::RemoveAtom(Atom& atom) {
    if (!HasAtom(atom)) return false;
    
    // Remove all bonds this atom is part of from molecule
    auto bnd_pred = [&atom](sBond bnd) {  // Predicate checks if atom in bnd
      auto s = bnd->GetAtoms();
      return !(&s.first == &atom || &s.second == &atom);
    };
    // Partition so can remove bonds from other atoms
    auto bnd_pos = std::partition(_bnds.begin(), _bnds.end(), bnd_pred);
    auto bnd_pos_erase = bnd_pos;
    for (auto it = _bnds.end(); bnd_pos != it; ++bnd_pos) {
      auto s = (*bnd_pos)->GetAtoms();
      if (&s.first == &atom) s.second.RemoveBond(*(*bnd_pos));
      else s.first.RemoveBond(*(*bnd_pos));
    }
    _bnds.erase(bnd_pos_erase, _bnds.end());
    
    // Remove all angles this atom is part of from molecule
    auto ang_pred = [&atom](sAngle ang) {
      auto s = ang->GetAtoms();
      return !(&s.first == &atom || &s.second == &atom || &s.third == &atom);
    };
    auto ang_pos = std::partition(_angs.begin(), _angs.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = _angs.end(); ang_pos != it; ++ang_pos) {
      auto s = (*ang_pos)->GetAtoms();
      if (&s.first == &atom) {
        s.second.RemoveAngle(*(*ang_pos));
        s.third.RemoveAngle(*(*ang_pos));
      } else if (&s.second == &atom) {
        s.first.RemoveAngle(*(*ang_pos));
        s.third.RemoveAngle(*(*ang_pos));
      } else {
        s.second.RemoveAngle(*(*ang_pos));
        s.third.RemoveAngle(*(*ang_pos));
      }
    }
    _angs.erase(ang_pos_erase, _angs.end());
    
    // Remove all dihedrals this atom is part of from molecule
    auto dhd_pred = [&atom](sDihedral dhd) {
      auto s = dhd->GetAtoms();
      return !(&s.first == &atom || &s.second == &atom
               || &s.third == &atom || &s.fourth == &atom);
    };
    auto dhd_pos = std::partition(_dhds.begin(), _dhds.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = _dhds.end(); dhd_pos != it; ++dhd_pos) {
      auto s = (*dhd_pos)->GetAtoms();
      if (&s.first == &atom) {
        s.second.RemoveDihedral(*(*dhd_pos));
        s.third.RemoveDihedral(*(*dhd_pos));
        s.fourth.RemoveDihedral(*(*dhd_pos));
      } else if (&s.second == &atom) {
        s.first.RemoveDihedral(*(*dhd_pos));
        s.third.RemoveDihedral(*(*dhd_pos));
        s.fourth.RemoveDihedral(*(*dhd_pos));
      } else if (&s.third == &atom) {
        s.second.RemoveDihedral(*(*dhd_pos));
        s.first.RemoveDihedral(*(*dhd_pos));
        s.fourth.RemoveDihedral(*(*dhd_pos));
      } else {
        s.second.RemoveDihedral(*(*dhd_pos));
        s.third.RemoveDihedral(*(*dhd_pos));
        s.first.RemoveDihedral(*(*dhd_pos));
      }
    }
    _dhds.erase(dhd_pos_erase, _dhds.end());
    
    // Remove the atom from the molecule
    graph::MGVertex v = _g->GetVertex(atom);
    _g->RemoveVertex(v);
    _atms.erase(std::find(_atms.begin(), _atms.end(), atom.shared_from_this()));
    ModificationMade();
    return true;
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule adding and removing atoms") {
    BuildBenzene();
    benzene.PerceiveAngles();
    benzene.PerceiveDihedrals();
    benzene.get_angs()[3]->SwapOrder(); // for better coverage of RemoveAtom
    
    Element test_element = GetPeriodicTable().GetElement("W");
    blankmol.Init();
    
    // New atom with no parameters
    Atom test_norm = blankmol.NewAtom();
    check_ne(Atom(), test_norm);
    check_eq(1, blankmol.NumAtoms());
    check_eq(blankmol.NumAtoms(), blankmol.GetGraph()->NumVertices());
    
    // New atom with element
    Atom test_elem = blankmol.NewAtom(test_element);
    check_ne(Atom(), test_elem);
    check_eq(2, blankmol.NumAtoms());
    check_eq(test_element, test_elem->GetElement());
    check_eq(blankmol.NumAtoms(), blankmol.GetGraph()->NumVertices());
    
    // New atom with name
    Atom test_name = blankmol.NewAtom("TestAtom");
    check_ne(Atom(), test_name);
    check_eq(3, blankmol.NumAtoms());
    check_eq("TestAtom", test_name->GetName());
    check_eq(blankmol.NumAtoms(), blankmol.GetGraph()->NumVertices());
    
    // New atom with name and element
    Atom test_name_elem = blankmol.NewAtom("TestAtomMore", test_element);
    check_ne(Atom(), test_name_elem);
    check_eq(4, blankmol.NumAtoms());
    check_eq("TestAtomMore", test_name_elem->GetName());
    check_eq(test_element, test_name_elem->GetElement());
    check_eq(blankmol.NumAtoms(), blankmol.GetGraph()->NumVertices());
    
    // Removing an atom
    check(blankmol.RemoveAtom(test_name_elem));
    check_eq(3, blankmol.NumAtoms());
    check_false(blankmol.HasAtom(test_name_elem));
    check_eq(Molecule(), test_name_elem->GetMolecule());
    check_eq(blankmol.NumAtoms(), blankmol.GetGraph()->NumVertices());
    
    // Fail cases
    check_false(blankmol.RemoveAtom(test_name_elem));
    check_false(blankmol.RemoveAtom(a_benzene.front()));
    check_false(blankmol.RemoveAtom(Atom()));
    
    // Removing an atom should remove all bonds, angles and dihedrals tied to it
    check(benzene.RemoveAtom(a_benzene.front()));
    check_eq(11, benzene.get_atms().size());
    check_eq(11, benzene.GetGraph()->NumVertices());
    check_eq(9, benzene.get_bnds().size());
    check_eq(9, benzene.GetGraph()->NumEdges());
    check_eq(11, benzene.get_angs().size());
    check_eq(12, benzene.get_dhds().size());
  }
 */
  
  Bond& Molecule::NewBond(Atom& a, Atom& b) {
    if (!HasAtom(a) || !HasAtom(b))
      throw std::out_of_range("Atoms not part of molecue");
    if (HasBond(a, b))
      throw std::out_of_range("Already has bond");
    sBond bn = std::make_shared<Bond>(a, b, *this);
    _bnds.emplace_back(bn);
    a.AddBond(*bn);
    b.AddBond(*bn);
    _g->AddEdge(*bn);
    ModificationMade();
    return *bn;
  }
  
  bool Molecule::RemoveBond(Bond& bond) {
    if (!HasBond(bond)) return false;
    Atom& a = bond.GetAtoms().first;
    Atom& b = bond.GetAtoms().second;
    // Remove the bond from each atom
    a.RemoveBond(bond);
    b.RemoveBond(bond);
    
    // Remove all angles this bond is part of from molecule
    auto ang_pred = [&a,&b] (sAngle ang) {
      auto s = ang->GetAtoms();
      return !((&a == &s.second && (&b == &s.first || &b == &s.third))
              || (&b == &s.second && (&a == &s.first || &a == &s.third)));
    };
    auto ang_pos = std::partition(_angs.begin(), _angs.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = _angs.end(); it != ang_pos; ++ang_pos) {
      auto s = (*ang_pos)->GetAtoms();
      s.first.RemoveAngle(*(*ang_pos));
      s.second.RemoveAngle(*(*ang_pos));
      s.third.RemoveAngle(*(*ang_pos));
    }
    _angs.erase(ang_pos_erase, _angs.end());
    
    // Remove all dihedrals this bond is part of from molecule
    auto dhd_pred = [&a,&b] (sDihedral dhd) {
      auto s = dhd->GetAtoms();
      return !((&a == &s.second && (&b == &s.first || &b == &s.third))
              || (&a == &s.third && (&b == &s.second || &b == &s.fourth))
              || (&b == &s.second && (&a == &s.first || &a == &s.third))
              || (&b == &s.third && (&a == &s.second || &a == &s.fourth)));
    };
    auto dhd_pos = std::partition(_dhds.begin(), _dhds.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = _dhds.end(); dhd_pos != it; ++dhd_pos) {
      auto s = (*dhd_pos)->GetAtoms();
      s.first.RemoveDihedral(*(*dhd_pos));
      s.second.RemoveDihedral(*(*dhd_pos));
      s.third.RemoveDihedral(*(*dhd_pos));
      s.fourth.RemoveDihedral(*(*dhd_pos));
    }
    _dhds.erase(dhd_pos_erase, _dhds.end());
    
    // Remove the bond from the molecule
    graph::MGEdge e = _g->GetEdge(bond);
    _g->RemoveEdge(e);
    _bnds.erase(std::find(_bnds.begin(), _bnds.end(), bond.shared_from_this()));
    ModificationMade();
    return true;
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule adding and removing bonds") {
    blankmol.Init();
    BuildBenzene();
    benzene.PerceiveAngles();
    benzene.PerceiveDihedrals();
    benzene.get_angs()[3]->SwapOrder(); // for better coverage of RemoveAtom
    benzene.get_dhds()[4]->SwapOrder(); // for better coverage of RemoveBond
    benzene.get_dhds()[0]->SwapOrder(); // for better coverage of RemoveBond
    
    Atom atm1 = blankmol.NewAtom();
    Atom atm2 = blankmol.NewAtom();
    Atom atm3 = blankmol.NewAtom();
    Atom atm4 = blankmol.NewAtom();
    
    // New bond
    Bond test_bnd =  blankmol.NewBond(atm1, atm2);
    check_ne(Bond(), test_bnd);
    check_eq(1, blankmol.NumBonds());
    check_eq(blankmol.NumBonds(), blankmol.GetGraph()->NumEdges());
    
    // New bond between existing should fail
    check_eq(Bond(), blankmol.NewBond(atm2, atm1));
    check_eq(1, blankmol.NumBonds());
    check_eq(blankmol.NumBonds(), blankmol.GetGraph()->NumEdges());
    
    // New bond with unowned atom should fail
    check_eq(Bond(), blankmol.NewBond(atm3, a_benzene.front()));
    check_eq(1, blankmol.NumBonds());
    check_eq(blankmol.NumBonds(), blankmol.GetGraph()->NumEdges());
    
    // New bond with null atom should fail
    check_eq(Bond(), blankmol.NewBond(atm4, Atom()));
    check_eq(1, blankmol.NumBonds());
    check_eq(blankmol.NumBonds(), blankmol.GetGraph()->NumEdges());
    Bond test_bnd2 = blankmol.NewBond(atm3, atm4);
    
    // Remove a bond
    check(blankmol.RemoveBond(test_bnd));
    check_eq(1, blankmol.NumBonds());
    check_eq(blankmol.NumBonds(), blankmol.GetGraph()->NumEdges());
    check_eq(Molecule(), test_bnd->GetMolecule());
    
    // Remove a bond between atoms
    check(blankmol.RemoveBond(atm4, atm3));
    check_eq(0, blankmol.NumBonds());
    check_eq(blankmol.NumBonds(), blankmol.GetGraph()->NumEdges());
    check_eq(Molecule(), test_bnd2->GetMolecule());
    
    // Fail cases
    check_false(blankmol.RemoveBond(test_bnd));
    check_false(blankmol.RemoveBond(Bond()));
    check_false(blankmol.RemoveBond(atm1, atm2));
    check_false(blankmol.RemoveBond(Atom(), atm3));
    
    // Removing a bond should remove all angles and dihedreals tied to it
    check(benzene.RemoveBond(b_benzene.front()));
    check_false(benzene.HasBond(a_benzene[0], a_benzene[1]));
    check_eq(12, benzene.get_atms().size());
    check_eq(11, benzene.get_bnds().size());
    check_eq(11, benzene.GetGraph()->NumEdges());
    check_eq(14, benzene.get_angs().size());
    check_eq(16, benzene.get_dhds().size());
  }
 */
  
  void Molecule::NewAngle(Atom &a, Atom &b, Atom &c) {
    // No need for logic checks as no ability for user to add angles
    _angs.emplace_back(std::make_shared<Angle>(a, b, c, *this));
    a.AddAngle(*_angs.back());
    b.AddAngle(*_angs.back());
    c.AddAngle(*_angs.back());
  }
  
  Dihedral& Molecule::NewDihedral(Atom &a, Atom &b, Atom &c, Atom &d) {
    if (!HasAtom(a) || !HasAtom(b) || !HasAtom(c) || !HasAtom(d))
      throw std::out_of_range("Atoms not part of molecule");
    if (HasDihedral(a, b, c, d))
      throw std::out_of_range("Already has dihedral");
    sDihedral dl = std::make_shared<Dihedral>(a, b, c, d, *this);
    _dhds.emplace_back(dl);
    a.AddDihedral(*dl);
    b.AddDihedral(*dl);
    c.AddDihedral(*dl);
    d.AddDihedral(*dl);
    ModificationMade();
    return *dl;
  }
  
/*  test_case_fixture(test::MoleculeTestFixture, "IXMolecule adding angles/dihedrals") {
    BuildRandomMolecule();
    check_ne(Angle(), randmol.NewAngle(a_randmol[0], a_randmol[1],
                                       a_randmol[2]));
    check_eq(1, randmol.get_angs().size());
    check_ne(Dihedral(), randmol.NewDihedral(a_randmol[0], a_randmol[2],
                                             a_randmol[3], a_randmol[4]));
    check_eq(1, randmol.get_dhds().size());
  }
 */
  
  sMolecule CreateMolecule() {
    sMolecule tmp = sMolecule(new Molecule());
    tmp->Init();
    return tmp;
  }
  
  test_suite_close();
}
