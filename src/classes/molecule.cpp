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
  test_suite_open("IXMolecule");
  
// ============================================================================
// == SERIALISATION ===========================================================
// ============================================================================
  
  template <typename Archive>
  void IXMolecule::Serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("atoms", _atms),
            INDIGOX_SERIAL_NVP("bonds", _bnds),
            INDIGOX_SERIAL_NVP("angles", _angs),
            INDIGOX_SERIAL_NVP("dihedrals", _dhds),
            INDIGOX_SERIAL_NVP("name", _name),
            INDIGOX_SERIAL_NVP("molecular_charge", _q),
            INDIGOX_SERIAL_NVP("molecular_graph", _g));
  }
  INDIGOX_SERIALISE(IXMolecule);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMolecule serialisation", T, ixmol_serial) {
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
  
// ============================================================================
// == CONSTRUCTION AND INITALISATION ==========================================
// ============================================================================
  
  IXMolecule::IXMolecule() : utils::IXCountableObject<IXMolecule>(), _name(""),
  _q(0), _emerge(0) {
    _emerge.set();
  }
  
  void IXMolecule::Init() {
    _g = graph::MolecularGraph(new graph::IXMolecularGraph(shared_from_this()));
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule construction") {
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
  
// ============================================================================
// == PROPERTY MODIFICATION ===================================================
// ============================================================================
  
  void IXMolecule::SetPropertyModified(Property prop) {
    EmergeSet mask(0);
    switch (prop) {
      case Property::ATOM_ELEMENTS:
        mask.set(static_cast<size_t>(Emergent::MOLECULAR_FORMULA));
        mask.set(static_cast<size_t>(Emergent::TOPOLOGICAL_BOFC));
        break;
      case Property::CONNECTIVITY:
        mask.set(static_cast<size_t>(Emergent::ANGLE_PERCEPTION));
        mask.set(static_cast<size_t>(Emergent::DIHEDRAL_PERCEPTION));
        mask.set(static_cast<size_t>(Emergent::TOPOLOGICAL_BOFC));
        break;
      case Property::ELECTRON_COUNT:
        mask.set(static_cast<size_t>(Emergent::TOPOLOGICAL_BOFC));
        break;
      case Property::NUM_PROPERTIES:
      default:
        break;
    }
    _emerge |= mask;
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule modify properties") {
    BuildRandomMolecule();
    size_t mol_form = static_cast<size_t>(EmergeProp::MOLECULAR_FORMULA);
    size_t top_bofc = static_cast<size_t>(EmergeProp::TOPOLOGICAL_BOFC);
    size_t ang_perc = static_cast<size_t>(EmergeProp::ANGLE_PERCEPTION);
    size_t dhd_perc = static_cast<size_t>(EmergeProp::DIHEDRAL_PERCEPTION);
    randmol.get_emerge().reset();
    
    subcase("property modifying") {
      // ATOM_ELEMENTS should set MOLECULAR_FORMULA and TOPOLOGICAL_BOFC
      randmol.SetPropertyModified(MolProperty::ATOM_ELEMENTS);
      check(randmol.get_emerge()[mol_form]);
      check(randmol.get_emerge()[top_bofc]);
      check_false(randmol.get_emerge()[ang_perc]);
      check_false(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // CONNECTIVITY should set TOPOLOGICAL_BOFC and perceptions
      randmol.SetPropertyModified(MolProperty::CONNECTIVITY);
      check_false(randmol.get_emerge()[mol_form]);
      check(randmol.get_emerge()[top_bofc]);
      check(randmol.get_emerge()[ang_perc]);
      check(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // ELECTRON_COUNT should set TOPOLOGICAL_BOFC
      randmol.SetPropertyModified(MolProperty::ELECTRON_COUNT);
      check_false(randmol.get_emerge()[mol_form]);
      check(randmol.get_emerge()[top_bofc]);
      check_false(randmol.get_emerge()[ang_perc]);
      check_false(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // NUM_PROPERTIES should set nothing
      randmol.SetPropertyModified(MolProperty::NUM_PROPERTIES);
      check_false(randmol.get_emerge()[mol_form]);
      check_false(randmol.get_emerge()[top_bofc]);
      check_false(randmol.get_emerge()[ang_perc]);
      check_false(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
    }
    
    subcase("molecular formula") {
      // Adding an atom should set
      Atom test_atm = randmol.NewAtom();
      check(randmol.get_emerge()[mol_form]);
      randmol.get_emerge().reset();
      
      // Removing an atom should set
      randmol.RemoveAtom(test_atm);
      check(randmol.get_emerge()[mol_form]);
      randmol.get_emerge().reset();
      
      // Changing an atom's element should set
      a_randmol.front()->SetElement("Zr");
      check(randmol.get_emerge()[mol_form]);
      randmol.get_emerge().reset();
    }
    
    subcase("topological bond order/formal charge") {
      // Adding an atom should set
      Atom test_atm = randmol.NewAtom();
      check(randmol.get_emerge()[top_bofc]);
      randmol.get_emerge().reset();
      
      // Adding a bond should set
      Bond test_bnd = randmol.NewBond(test_atm, a_randmol.back());
      check(randmol.get_emerge()[top_bofc]);
      randmol.get_emerge().reset();
      
      // Removing a bond should set
      check(randmol.RemoveBond(b_randmol.front()));
      check(randmol.get_emerge()[top_bofc]);
      randmol.get_emerge().reset();
      
      // Removing a bond between atoms should set
      check(randmol.RemoveBond(a_randmol.back(), test_atm));
      check(randmol.get_emerge()[top_bofc]);
      randmol.get_emerge().reset();
      
      // Removing an atom should set
      check(randmol.RemoveAtom(test_atm));
      check(randmol.get_emerge()[top_bofc]);
      randmol.get_emerge().reset();
      
      // Changing an element should set
      a_randmol.front()->SetElement("W");
      check(randmol.get_emerge()[top_bofc]);
      randmol.get_emerge().reset();
      
      // Changing the molecular charge should set
      randmol.SetMolecularCharge(-23);
      check(randmol.get_emerge()[top_bofc]);
      randmol.get_emerge().reset();
    }
    
    subcase("angle perception") {
      // Adding an atom shouldn't set
      Atom test_atm = randmol.NewAtom();
      check_false(randmol.get_emerge()[ang_perc]);
      randmol.get_emerge().reset();
      
      // Adding a bond should set
      Bond test_bnd = randmol.NewBond(test_atm, a_randmol.back());
      check(randmol.get_emerge()[ang_perc]);
      randmol.get_emerge().reset();
      
      // Removing a bond should set
      check(randmol.RemoveBond(b_randmol.front()));
      check(randmol.get_emerge()[ang_perc]);
      randmol.get_emerge().reset();
      
      // Removing a bond between atoms should set
      check(randmol.RemoveBond(test_atm, a_randmol.back()));
      check(randmol.get_emerge()[ang_perc]);
      randmol.get_emerge().reset();
      
      // Removing an atom without bonds should not set
      check(randmol.RemoveAtom(test_atm));
      check_false(randmol.get_emerge()[ang_perc]);
      randmol.get_emerge().reset();
      
      // Removing an atom with bonds should set
      for (Atom atm : a_randmol) {
        if (atm->NumBonds()) {
          check(randmol.RemoveAtom(atm));
          check(randmol.get_emerge()[ang_perc]);
          randmol.get_emerge().reset();
          break;
        }
      }
    }
    
    subcase("dihedral perception") {
      // Adding an atom shouldn't set
      Atom test_atm = randmol.NewAtom();
      check_false(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // Adding a bond should set
      Bond test_bnd = randmol.NewBond(test_atm, a_randmol.back());
      check(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // Removing a bond should set
      check(randmol.RemoveBond(b_randmol.front()));
      check(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // Removing a bond between atoms should set
      check(randmol.RemoveBond(test_atm, a_randmol.back()));
      check(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // Removing an atom without bonds should not set
      check(randmol.RemoveAtom(test_atm));
      check_false(randmol.get_emerge()[dhd_perc]);
      randmol.get_emerge().reset();
      
      // Removing an atom with bonds should set
      for (Atom atm : a_randmol) {
        if (atm->NumBonds()) {
          check(randmol.RemoveAtom(atm));
          check(randmol.get_emerge()[dhd_perc]);
          randmol.get_emerge().reset();
          break;
        }
      }
    }
  }
  
// ============================================================================
// == FINDING ITEMS ===========================================================
// ============================================================================
  
  Bond IXMolecule::_FindBond(const Atom &a, const Atom &b) const {
    auto pred = [a,b](Bond bnd) {
      std::pair<Atom, Atom> atms = bnd->GetAtoms();
      return ((atms.first == a && atms.second == b)
              || (atms.first == b && atms.second == a));
    };
    if (!a || !b) return Bond();
    MolBondIter pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return pos == _bnds.end() ? Bond() : *pos;
  }
  
  Angle IXMolecule::_FindAngle(const Atom &a, const Atom &b, const Atom &c) const {
    stdx::triple<Atom, Atom, Atom> a2c = stdx::make_triple(a, b, c);
    stdx::triple<Atom, Atom, Atom> c2a = stdx::make_triple(c, b, a);
    auto pred = [&a2c, &c2a] (Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return (atms == a2c || atms == c2a);
    };
    MolAngleIter pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return pos == _angs.end() ? Angle() : *pos;
  }
  
  Dihedral IXMolecule::_FindDihedral(const Atom &a, const Atom &b,
                                     const Atom &c, const Atom &d) const {
    stdx::quad<Atom, Atom, Atom, Atom> a2d = stdx::make_quad(a, b, c, d);
    stdx::quad<Atom, Atom, Atom, Atom> d2a = stdx::make_quad(d, c, b, a);
    auto pred = [&a2d, &d2a] (Dihedral dhd) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = dhd->GetAtoms();
      return (atms == a2d || atms == d2a);
    };
    MolDihedralIter pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    return pos == _dhds.end() ? Dihedral() : *pos;
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule find items") {
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
  
  bool IXMolecule::HasBond(const Atom& a, const Atom& b) const {
    if (!HasAtom(a) || !HasAtom(b)) return false;
    return bool(_FindBond(a, b));
  }
  
  bool IXMolecule::HasAngle(const Atom &a, const Atom &b, const Atom &c) {
    if (!HasAtom(a) || !HasAtom(b) || !HasAtom(c)) return false;
    PerceiveAngles();
    return bool(_FindAngle(a, b, c));
  }
  
  bool IXMolecule::HasDihedral(const Atom &a, const Atom &b,
                               const Atom &c, const Atom &d) {
    if (!HasAtom(a) || !HasAtom(b) || !HasAtom(c) || !HasAtom(d)) return false;
    PerceiveDihedrals();
    return bool(_FindDihedral(a, b, c, d));
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule has items") {
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
  
// ============================================================================
// == GETTING AND SETTING =====================================================
// ============================================================================
  
  Atom IXMolecule::GetAtomTag(uint32_t tag) const {
    auto pred = [tag](Atom a) { return a->GetTag() == tag; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? Atom() : *pos;
  }
  
  Atom IXMolecule::GetAtomID(uint32_t id) const {
    auto pred = [id](Atom a) { return a->GetUniqueID() == id; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? Atom() : *pos;
  }
  
  Bond IXMolecule::GetBond(const Atom& a, const Atom& b) const {
    return _FindBond(a, b);
  }
  
  Bond IXMolecule::GetBondTag(uint32_t tag) const {
    auto pred = [tag](Bond b) { return b->GetTag() == tag; };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return (pos == _bnds.end()) ? Bond() : *pos;
  }
  
  Bond IXMolecule::GetBondID(uint32_t id) const {
    auto pred = [id](Bond b) { return b->GetUniqueID() == id; };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return (pos == _bnds.end()) ? Bond() : *pos;
  }
  
  Angle IXMolecule::GetAngle(const Atom& a, const Atom&b, const Atom& c) {
    PerceiveAngles();
    return _FindAngle(a, b, c);
  }
  
  Angle IXMolecule::GetAngleTag(uint32_t tag) const {
    auto pred = [tag](Angle ang) { return ang->GetTag() == tag; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return (pos == _angs.end()) ? Angle() : *pos;
  }
  
  Dihedral IXMolecule::GetDihedral(const Atom &a, const Atom &b,
                                   const Atom &c, const Atom &d) {
    PerceiveDihedrals();
    return _FindDihedral(a, b, c, d);
  }
  
  Dihedral IXMolecule::GetDihedralTag(uint32_t tag) const {
    auto pred = [tag](Dihedral dhd) { return dhd->GetTag() == tag; };
    auto pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    return (pos == _dhds.end()) ? Dihedral() : *pos;
  }
  
  Angle IXMolecule::GetAngleID(uint32_t id) const {
    auto pred = [id](Angle ang) { return ang->GetUniqueID() == id; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return (pos == _angs.end()) ? Angle() : *pos;
  }
  
  Dihedral IXMolecule::GetDihedralID(uint32_t id) const {
    auto pred = [id](Dihedral dhd) { return dhd->GetUniqueID() == id; };
    auto pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    return (pos == _dhds.end()) ? Dihedral() : *pos;
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule getting and setting") {
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
  
  std::string IXMolecule::GetFormula() {
    if (_emerge[static_cast<size_t>(Emergent::MOLECULAR_FORMULA)]) {
      std::map<std::string, size_t> e_count;
      for (Atom atm : _atms) e_count[atm->GetElement()->GetSymbol()]++;
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
      _formula_cache = ss.str();
      _emerge.reset(static_cast<size_t>(Emergent::MOLECULAR_FORMULA));
    }
    return _formula_cache;
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule getting formula") {
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
  
// ============================================================================
// == PERCEPTION ==============================================================
// ============================================================================
  
  size_t IXMolecule::PerceiveAngles() {
    using namespace indigox::graph;
    if (!_emerge[static_cast<size_t>(Emergent::ANGLE_PERCEPTION)]) return 0;
    
    // Expected number of angles
    auto sum = [&](size_t current, Atom v) -> size_t {
      size_t degree = v->NumBonds();
      if (degree < 2) return current;
      return current + degree * (degree - 1) / 2;
    };
    size_t count = std::accumulate(_atms.begin(), _atms.end(), 0, sum);
    _angs.reserve(count);
    count = 0;
    
    // Adding new angles
    for (MolAtomIter b = _atms.begin(), e = _atms.end(); b != e; ++b) {
      if ((*b)->NumBonds() < 2) continue;
      auto nbrs_it = _g->GetNeighbours(_g->GetVertex(*b));
      std::vector<MGVertex> nbrs(nbrs_it.first, nbrs_it.second);
      for (size_t i = 0; i < nbrs.size() - 1; ++i) {
        Atom a = nbrs[i]->GetAtom();
        for (size_t j = i + 1; j < nbrs.size(); ++j) {
          Atom c = nbrs[j]->GetAtom();
          if (_FindAngle(a, *b, c)) continue;
          NewAngle(a, *b, c);
          ++count;
        }
      }
    }
    _emerge.reset(static_cast<size_t>(Emergent::ANGLE_PERCEPTION));
    return count;
  }
  
  size_t IXMolecule::PerceiveDihedrals() {
    using namespace indigox::graph;
    if (!_emerge[static_cast<size_t>(Emergent::DIHEDRAL_PERCEPTION)]) return 0;
    
    // Expected number of dihedrals
    auto sum = [&](size_t current, Bond b) -> size_t {
      size_t b_degree = b->GetSourceAtom()->NumBonds();
      if (b_degree < 2) return current;
      size_t c_degree = b->GetTargetAtom()->NumBonds();
      if (c_degree < 2) return current;
      return current + (b_degree - 1) * (c_degree - 1);
    };
    size_t count = std::accumulate(_bnds.begin(), _bnds.end(), 0, sum);
    _dhds.reserve(count);
    count = 0;
    
    // Adding new dihedrals
    for (MolBondIter b = _bnds.begin(), e = _bnds.end(); b != e; ++b) {
      Atom B, C;
      std::tie(B,C) = (*b)->GetAtoms();
      if (B->NumBonds() < 2 || C->NumBonds() < 2) continue;
      auto nbrs_it = _g->GetNeighbours(_g->GetVertex(B));
      std::vector<MGVertex> As(nbrs_it.first, nbrs_it.second);
      nbrs_it = _g->GetNeighbours(_g->GetVertex(C));
      std::vector<MGVertex> Ds(nbrs_it.first, nbrs_it.second);
      for (size_t i = 0; i < As.size(); ++i) {
        Atom A = As[i]->GetAtom();
        if (A == C) continue;
        for (size_t j = 0; j < Ds.size(); ++j) {
          Atom D = Ds[j]->GetAtom();
          if (D == B) continue;
          if (_FindDihedral(A, B, C, D)) continue;
          NewDihedral(A, B, C, D);
          ++count;
        }
      }
    }
    _emerge.reset(static_cast<size_t>(Emergent::DIHEDRAL_PERCEPTION));
    return count;
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule perception") {
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
 
// ============================================================================
// == STRUCTURE MODIFICATION ==================================================
// ============================================================================
  
  Atom IXMolecule::NewAtom() {
    Atom atom = std::make_shared<IXAtom>(shared_from_this());
    _atms.emplace_back(atom);
    _g->AddVertex(atom);
    SetPropertyModified(Property::ATOM_ELEMENTS);
    return atom;
  }
  
  Atom IXMolecule::NewAtom(Element element) {
    Atom atom = NewAtom();
    atom->SetElement(element);
    return atom;
  }
  
  Atom IXMolecule::NewAtom(std::string name) {
    Atom atom = NewAtom();
    atom->SetName(name);
    return atom;
  }
  
  Atom IXMolecule::NewAtom(std::string name, Element element) {
    Atom atom = NewAtom();
    atom->SetName(name);
    atom->SetElement(element);
    return atom;
  }
  
  bool IXMolecule::RemoveAtom(Atom atom) {
    if (!atom || !HasAtom(atom)) return false;
    
    // Remove all bonds this atom is part of from molecule
    if (atom->NumBonds()) SetPropertyModified(Property::CONNECTIVITY);
    auto bnd_pred = [atom](Bond bnd) {  // Predicate checks if atom in bnd
      std::pair<Atom, Atom> atms = bnd->GetAtoms();
      return !(atms.first == atom || atms.second == atom);
    };
    // Partition so can remove bonds from other atoms
    auto bnd_pos = std::partition(_bnds.begin(), _bnds.end(), bnd_pred);
    auto bnd_pos_erase = bnd_pos;
    for (auto it = _bnds.end(); bnd_pos != it; ++bnd_pos) {
      std::pair<Atom, Atom> atms = (*bnd_pos)->GetAtoms();
      if (atms.first == atom) atms.second->RemoveBond(*bnd_pos);
      else atms.first->RemoveBond(*bnd_pos);
    }
    _bnds.erase(bnd_pos_erase, _bnds.end());
    
    // Remove all angles this atom is part of from molecule
    auto ang_pred = [atom](Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return !(atms.first == atom || atms.second == atom || atms.third == atom);
    };
    auto ang_pos = std::partition(_angs.begin(), _angs.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = _angs.end(); ang_pos != it; ++ang_pos) {
      stdx::triple<Atom, Atom, Atom> atms = (*ang_pos)->GetAtoms();
      if (atms.first == atom) {
        atms.second->RemoveAngle(*ang_pos);
        atms.third->RemoveAngle(*ang_pos);
      } else if (atms.second == atom) {
        atms.first->RemoveAngle(*ang_pos);
        atms.third->RemoveAngle(*ang_pos);
      } else {
        atms.second->RemoveAngle(*ang_pos);
        atms.third->RemoveAngle(*ang_pos);
      }
    }
    _angs.erase(ang_pos_erase, _angs.end());
    
    // Remove all dihedrals this atom is part of from molecule
    auto dhd_pred = [atom](Dihedral dhd) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = dhd->GetAtoms();
      return !(atms.first == atom || atms.second == atom
               || atms.third == atom || atms.fourth == atom);
    };
    auto dhd_pos = std::partition(_dhds.begin(), _dhds.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = _dhds.end(); dhd_pos != it; ++dhd_pos) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = (*dhd_pos)->GetAtoms();
      if (atms.first == atom) {
        atms.second->RemoveDihedral(*dhd_pos);
        atms.third->RemoveDihedral(*dhd_pos);
        atms.fourth->RemoveDihedral(*dhd_pos);
      } else if (atms.second == atom) {
        atms.first->RemoveDihedral(*dhd_pos);
        atms.third->RemoveDihedral(*dhd_pos);
        atms.fourth->RemoveDihedral(*dhd_pos);
      } else if (atms.third == atom) {
        atms.second->RemoveDihedral(*dhd_pos);
        atms.first->RemoveDihedral(*dhd_pos);
        atms.fourth->RemoveDihedral(*dhd_pos);
      } else {
        atms.second->RemoveDihedral(*dhd_pos);
        atms.third->RemoveDihedral(*dhd_pos);
        atms.first->RemoveDihedral(*dhd_pos);
      }
    }
    _dhds.erase(dhd_pos_erase, _dhds.end());
    
    // Remove the atom from the molecule
    _atms.erase(std::find(_atms.begin(), _atms.end(), atom));
    _g->RemoveVertex(_g->GetVertex(atom));
    SetPropertyModified(Property::ATOM_ELEMENTS);
    
    // Clear the atom to invalidate it
    atom->Clear();
    return true;
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule adding and removing atoms") {
    BuildBenzene();
    benzene.PerceiveAngles();
    benzene.PerceiveDihedrals();
    benzene.get_angs()[3]->SwapOrder(); // for better coverage of RemoveAtom
    
    Element test_element = GetPeriodicTable()->GetElement("W");
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
  
  Bond IXMolecule::NewBond(Atom a, Atom b) {
    if (!a || !b || !HasAtom(a) || !HasAtom(b)) return Bond();
    Bond bond = GetBond(a, b);
    if (bond) return Bond();
    bond = std::make_shared<IXBond>(a, b, shared_from_this());
    _bnds.push_back(bond);
    a->AddBond(bond);
    b->AddBond(bond);
    _g->AddEdge(bond);
    SetPropertyModified(Property::CONNECTIVITY);
    return  bond;
  }
  
  bool IXMolecule::RemoveBond(Bond bond) {
    if (!bond || !HasBond(bond)) return false;
    Atom a, b;
    std::tie(a, b) = bond->GetAtoms();
    // Remove the bond from each atom
    a->RemoveBond(bond);
    b->RemoveBond(bond);
    
    // Remove all angles this bond is part of from molecule
    auto ang_pred = [a,b] (Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return !((a == atms.second && (b == atms.first || b == atms.third))
              || (b == atms.second && (a == atms.first || a == atms.third)));
    };
    auto ang_pos = std::partition(_angs.begin(), _angs.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = _angs.end(); it != ang_pos; ++ang_pos) {
      stdx::triple<Atom, Atom, Atom> atms = (*ang_pos)->GetAtoms();
      atms.first->RemoveAngle(*ang_pos);
      atms.second->RemoveAngle(*ang_pos);
      atms.third->RemoveAngle(*ang_pos);
    }
    _angs.erase(ang_pos_erase, _angs.end());
    
    // Remove all dihedrals this bond is part of from molecule
    auto dhd_pred = [a,b] (Dihedral dhd) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = dhd->GetAtoms();
      return !((a == atms.second && (b == atms.first || b == atms.third))
              || (a == atms.third && (b == atms.second || b == atms.fourth))
              || (b == atms.second && (a == atms.first || a == atms.third))
              || (b == atms.third && (a == atms.second || a == atms.fourth)));
    };
    auto dhd_pos = std::partition(_dhds.begin(), _dhds.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = _dhds.end(); dhd_pos != it; ++dhd_pos) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = (*dhd_pos)->GetAtoms();
      atms.first->RemoveDihedral(*dhd_pos);
      atms.second->RemoveDihedral(*dhd_pos);
      atms.third->RemoveDihedral(*dhd_pos);
      atms.fourth->RemoveDihedral(*dhd_pos);
    }
    _dhds.erase(dhd_pos_erase, _dhds.end());
    
    // Remove the bond from the molecule
    _bnds.erase(std::find(_bnds.begin(), _bnds.end(), bond));
    _g->RemoveEdge(_g->GetEdge(bond));
    SetPropertyModified(Property::CONNECTIVITY);
    
    // Clear the bond to invalidate it
    bond->Clear();
    return true;
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule adding and removing bonds") {
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
  
  Angle IXMolecule::NewAngle(const Atom &a, const Atom &b, const Atom &c) {
    // No need for logic checks as no ability for user to add angles
    _angs.emplace_back(std::make_shared<IXAngle>(a, b, c, shared_from_this()));
    a->AddAngle(_angs.back());
    b->AddAngle(_angs.back());
    c->AddAngle(_angs.back());
    return _angs.back();
  }
  
  Dihedral IXMolecule::NewDihedral(const Atom &a, const Atom &b,
                                   const Atom &c, const Atom &d) {
    // No need for logic checks as no ability for user to add dihedrals
    _dhds.emplace_back(std::make_shared<IXDihedral>(a, b, c, d, shared_from_this()));
    a->AddDihedral(_dhds.back());
    b->AddDihedral(_dhds.back());
    c->AddDihedral(_dhds.back());
    d->AddDihedral(_dhds.back());
    return _dhds.back();
  }
  
  test_case_fixture(test::MoleculeTestFixture, "IXMolecule adding angles/dihedrals") {
    BuildRandomMolecule();
    check_ne(Angle(), randmol.NewAngle(a_randmol[0], a_randmol[1],
                                       a_randmol[2]));
    check_eq(1, randmol.get_angs().size());
    check_ne(Dihedral(), randmol.NewDihedral(a_randmol[0], a_randmol[2],
                                             a_randmol[3], a_randmol[4]));
    check_eq(1, randmol.get_dhds().size());
  }
  
  test_suite_close();
}
