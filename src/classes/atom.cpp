#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/atom_test.hpp>
#include <indigox/test/angle_test.hpp>
#include <indigox/test/bond_test.hpp>
#include <indigox/test/dihedral_test.hpp>
#include <indigox/test/forcefield_test.hpp>

namespace indigox {
  
  test::AtomTestFixture::AtomTestFixture()
  : mol(CreateMolecule()), atm(mol), fftype(CreateGenericTestFFAtom().imp) { }
  
  test_suite_open("IXAtom");
  
  template<typename Archive>
  void IXAtom::save(Archive &archive, const uint32_t) const {
    string_ element = GetElement()->GetName();
    archive(INDIGOX_SERIAL_NVP("molecule", _mol),
            INDIGOX_SERIAL_NVP("element", element),
            INDIGOX_SERIAL_NVP("formal_charge", _fc),
            INDIGOX_SERIAL_NVP("tag", _tag),
            INDIGOX_SERIAL_NVP("implicit_h_count", _implicitH),
            INDIGOX_SERIAL_NVP("name", _name),
            INDIGOX_SERIAL_NVP("position_x", _pos[0]),
            INDIGOX_SERIAL_NVP("position_y", _pos[1]),
            INDIGOX_SERIAL_NVP("position_z", _pos[2]),
            INDIGOX_SERIAL_NVP("partial_charge", _partial),
            INDIGOX_SERIAL_NVP("stereochemistry", _stereo),
            INDIGOX_SERIAL_NVP("is_aromatic", _aromatic),
            INDIGOX_SERIAL_NVP("type", _type),
            INDIGOX_SERIAL_NVP("bonds", _bnds),
            INDIGOX_SERIAL_NVP("angles", _angs),
            INDIGOX_SERIAL_NVP("dihedrals", _dhds)
            );
  }
  
  template <typename Archive>
  void IXAtom::load_and_construct(Archive &archive,
                                  cereal::construct<IXAtom> &construct,
                                  const uint32_t) {
    Molecule mol;
    string_ element;
    archive(INDIGOX_SERIAL_NVP("molecule", mol));
    construct(mol);
    archive(INDIGOX_SERIAL_NVP("element", element),
            INDIGOX_SERIAL_NVP("formal_charge", construct->_fc),
            INDIGOX_SERIAL_NVP("tag", construct->_tag),
            INDIGOX_SERIAL_NVP("implicit_h_count", construct->_implicitH),
            INDIGOX_SERIAL_NVP("name", construct->_name),
            INDIGOX_SERIAL_NVP("position_x", construct->_pos(0)),
            INDIGOX_SERIAL_NVP("position_y", construct->_pos(1)),
            INDIGOX_SERIAL_NVP("position_z", construct->_pos(2)),
            INDIGOX_SERIAL_NVP("partial_charge", construct->_partial),
            INDIGOX_SERIAL_NVP("stereochemistry", construct->_stereo),
            INDIGOX_SERIAL_NVP("is_aromatic", construct->_aromatic),
            INDIGOX_SERIAL_NVP("type", construct->_type),
            INDIGOX_SERIAL_NVP("bonds", construct->_bnds),
            INDIGOX_SERIAL_NVP("angles", construct->_angs),
            INDIGOX_SERIAL_NVP("dihedrals", construct->_dhds));
    construct->SetElement(element);
  }
  INDIGOX_SERIALISE_SPLIT(IXAtom);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXAtom serialisation", T, ixatom_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
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
      check_nothrow(oar(saved, bond, angle, dihedral, fixture.mol, fixture.fftype));
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
    check_eq(saved->GetVector(), loaded->GetVector());
    check_eq(saved->GetPartialCharge(), loaded->GetPartialCharge());
    check_eq(saved->GetStereochemistry(), loaded->GetStereochemistry());
    check_eq(saved->GetAromaticity(), loaded->GetAromaticity());
    check_eq(bond->GetTag(), loaded->GetBondIters().first->lock()->GetTag());
    check_eq(angle->GetTag(), loaded->GetAngleIters().first->lock()->GetTag());
    check_eq(dihedral->GetTag(), loaded->GetDihedralIters().first->lock()->GetTag());
    check_eq(bond_loaded, fixture.atm.get_bnds().front().lock());
    check_eq(angle_loaded, fixture.atm.get_angs().front().lock());
    check_eq(dihedral_loaded, fixture.atm.get_dhds().front().lock());
    check_eq(fftype_loaded, loaded->GetType());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixatom_serial, ixserial<IXAtom>);
  
  inline void __set_property_modified(_Molecule mol, MolProperty p) {
    if (!mol.expired()) mol.lock()->SetPropertyModified(p);
  }
  
  IXAtom::IXAtom(Molecule m) : utils::IXCountableObject<IXAtom>(), _mol(m),
  _elem(Element()), _fc(0), _tag(0), _implicitH(0), _name(), _pos(0.0,0.0,0.0),
  _partial(0.0), _stereo(Stereo::UNDEFINED), _aromatic(false) { }
  
  test_case_fixture(test::AtomTestFixture, "IXAtom construction") {
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
  
  void IXAtom::SetElement(Element e) {
    if (e != GetElement()) {
      _elem = e;
      __set_property_modified(_mol, MolProperty::ATOM_ELEMENTS);
    }
  }
  
  size_ IXAtom::GetIndex() const {
    Molecule mol = _mol.lock();
    if (!mol) return GetTag();
    auto be = mol->GetAtoms();
    auto pos = std::find(be.first, be.second, shared_from_this());
    if (pos == be.second) return GetTag();
    return std::distance(be.first, pos);
  }
  
  test_case_fixture(test::AtomTestFixture, "IXAtom getting and setting") {
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
    check_eq(Vec3(4.5, -7.9002, 0.00004), atm.GetVector());
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
  
  test_case_fixture(test::AtomTestFixture, "IXAtom adding and removing") {
    /*  Remember, normal use is only through IXMolecule, so some assumptions
     *  are made. These are that will never attempt to add the same item twice,
     *  and will never attempt to remove an item which has not been added. */
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<size_> distribution(7,17);
    size_ num = distribution(generator);
    
    subcase("implicit hydrogen") {
      size_ count = 0;
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
      for (size_ i = 1; i <= num; ++i) {
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
      for (size_ i = 1; i <= num; ++i) {
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
      for (size_ i = 1; i <= num; ++i) {
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
      for (_Dihedral dhd : atm.get_dhds()) saved_dihedrals.emplace_back(dhd.lock());
      check_eq(dihedrals.size(), saved_dihedrals.size());
      check_eq(dihedrals, saved_dihedrals);
      
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
  
  string_ IXAtom::ToString() {
    std::stringstream ss;
    ss << "Atom(" << GetIndex() << ", " << GetElement()->GetSymbol() << ")";
    return ss.str();
  }
  
  test_case_fixture(test::AtomTestFixture, "IXAtom printing methods") {
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
  
  void IXAtom::Clear() {
    _mol.reset();
    _elem.reset();
    _fc = 0;
    _tag = 0;
    _implicitH = 0;
    _name = "";
    _pos << 0.0, 0.0, 0.0;
    _partial = 0.0;
    _stereo = Stereo::UNDEFINED;
    _aromatic = false;
    _type.reset();
    _bnds.clear();
    _angs.clear();
    _dhds.clear();
  }
  
  test_case_fixture(test::AtomTestFixture, "IXAtom clearing methods") {
    atm.SetElement("Cr");
    atm.SetFormalCharge(-3);
    atm.SetTag(12);
    atm.SetImplicitCount(4);
    atm.SetName("Testable");
    atm.SetPosition(1.2, 3.4, 5.6);
    atm.SetPartialCharge(-0.009);
    atm.SetStereochemistry(AtomStereo::R);
    atm.SetAromaticity(true);
    atm.SetType(fftype);
    atm.AddBond(test::CreateGenericTestBond().imp);
    atm.AddAngle(test::CreateGenericTestAngle().imp);
    atm.AddDihedral(test::CreateGenericTestDihedral().imp);
    
    // Pre checks
    check_eq(mol, atm.get_mol().lock());
    check_eq(GetPeriodicTable()->GetElement("Cr"), atm.get_elem().lock());
    check_eq(12, atm.get_tag());
    check_eq(4, atm.get_implicitH());
    check_eq("Testable", atm.get_name());
    check_eq(Vec3(1.2, 3.4, 5.6), atm.get_pos());
    check_eq(approximately(-0.009), atm.get_partial());
    check_eq(AtomStereo::R, atm.get_stereo());
    check_eq(true, atm.get_aromatic());
    check_eq(fftype, atm.get_type());
    check_eq(1, atm.get_bnds().size());
    check_eq(1, atm.get_angs().size());
    check_eq(1, atm.get_dhds().size());
    
    atm.Clear();
    // Post checks
    check_ne(mol, atm.get_mol().lock());
    check_ne(GetPeriodicTable()->GetElement("Cr"), atm.get_elem().lock());
    check_ne(12, atm.get_tag());
    check_ne(4, atm.get_implicitH());
    check_ne("Testable", atm.get_name());
    check_ne(Vec3(1.2, 3.4, 5.6), atm.get_pos());
    check_ne(approximately(-0.009), atm.get_partial());
    check_ne(AtomStereo::R, atm.get_stereo());
    check_ne(true, atm.get_aromatic());
    check_ne(fftype, atm.get_type());
    check_ne(1, atm.get_bnds().size());
    check_ne(1, atm.get_angs().size());
    check_ne(1, atm.get_dhds().size());
    check_eq(Molecule(), atm.get_mol().lock());
    check_eq(Element(), atm.get_elem().lock());
    check_eq(0, atm.get_tag());
    check_eq(0, atm.get_implicitH());
    check_eq("", atm.get_name());
    check_eq(Vec3(), atm.get_pos());
    check_eq(0.0, atm.get_partial());
    check_eq(AtomStereo::UNDEFINED, atm.get_stereo());
    check_eq(false, atm.get_aromatic());
    check_eq(FFAtom(), atm.get_type());
    check_eq(0, atm.get_bnds().size());
    check_eq(0, atm.get_angs().size());
    check_eq(0, atm.get_dhds().size());
  }
  
  test_suite_close();
}
