#include <sstream>
#include <stdexcept>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/bond_test.hpp>

namespace indigox {
  using namespace indigox::utils; // for IXCountableObject and WeakContains
  
  test_suite_open("IXBond");
  
  // Serialisation of IXBond
  template <typename Archive>
  void IXBond::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("molecule", _mol));
    std::vector<Atom> atoms;
    for (_Atom at : _atms) atoms.emplace_back(at.lock());
    archive(INDIGOX_SERIAL_NVP("atoms", atoms));
    archive(INDIGOX_SERIAL_NVP("tag", _tag),
            INDIGOX_SERIAL_NVP("bond_order", _order),
            INDIGOX_SERIAL_NVP("is_aromatic", _aromatic),
            INDIGOX_SERIAL_NVP("stereochemistry", _stereo)
            );
  }
  
  template <typename Archive>
  void IXBond::load_and_construct(Archive &archive,
                                  cereal::construct<IXBond> &construct,
                                  const uint32_t) {
    Molecule m;
    archive(INDIGOX_SERIAL_NVP("molecule", m));
    std::vector<Atom> atoms;
    archive(INDIGOX_SERIAL_NVP("atoms", atoms));
    construct(atoms[0], atoms[1], m);
    archive(INDIGOX_SERIAL_NVP("tag", construct->_tag),
            INDIGOX_SERIAL_NVP("bond_order", construct->_order),
            INDIGOX_SERIAL_NVP("is_aromatic", construct->_aromatic),
            INDIGOX_SERIAL_NVP("stereochemistry", construct->_stereo)
            );
  }
  INDIGOX_SERIALISE_SPLIT(IXBond);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXBond serialisation", T, id) {
    using In = typename T::t1;
    using Out = typename T::t2;
    test::BondTestFixture fixture;
    Bond saved = fixture.bnd.imp;
    
    saved->SetTag(12);
    saved->SetOrder(BondOrder::TRIPLE);
    saved->SetAromaticity(true);
    saved->SetStereochemistry(BondStereo::Z);
    
    std::ostringstream os;
    {
      Out oar(os);
      // atoms and molecule saved so they can be kept alive on load
      check_nothrow(oar(saved, saved->GetSourceAtom(), saved->GetTargetAtom(),
                        saved->GetMolecule()));
    }
    
    Bond loaded;
    Atom source, target;
    Molecule mol;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, source, target, mol));
    }
    fixture.bnd.imp = loaded;
    check(!fixture.bnd.get_atms()[0].expired());
    check(!fixture.bnd.get_atms()[1].expired());
    check(!fixture.bnd.get_mol().expired());
    check_eq(saved->GetTag(), loaded->GetTag());
    check_eq(saved->GetOrder(), loaded->GetOrder());
    check_eq(saved->GetAromaticity(), loaded->GetAromaticity());
    check_eq(saved->GetStereochemistry(), loaded->GetStereochemistry());
    check_eq(source->GetTag(), saved->GetSourceAtom()->GetTag());
    check_eq(target->GetTag(), saved->GetTargetAtom()->GetTag());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(id, ixserial<IXBond>);
  
  IXBond::IXBond(Atom a, Atom b, Molecule m) : IXCountableObject<IXBond>(),
  _mol(m), _tag(0), _order(Order::UNDEFINED), _aromatic(false),
  _stereo(Stereo::UNDEFINED), _atms({{a,b}}) { }
  
  test_case_fixture(test::BondTestFixture, "IXBond construction") {
    check_nothrow(test::TestBond(a,b,mol));
    check_eq(mol, bnd.get_mol().lock());
    check_eq(0, bnd.get_tag());
    check_eq(BondOrder::UNDEFINED, bnd.get_order());
    check_eq(false, bnd.get_aromatic());
    check_eq(BondStereo::UNDEFINED, bnd.get_stereo());
    check_eq(a, bnd.get_atms()[0].lock());
    check_eq(b, bnd.get_atms()[1].lock());
    
    // Check unique IDs correctly update
    test::TestBond bnd1(a,b,mol);
    test::TestBond bnd2(a,b,mol);
    check_ne(bnd1.GetUniqueID(), bnd2.GetUniqueID());
    check_eq(bnd1.GetUniqueID() + 1, bnd2.GetUniqueID());
  }
  
  test_case_fixture(test::BondTestFixture, "IXBond getting and setting") {
    // Check no throwing
    check_nothrow(bnd.GetTag());
    check_nothrow(bnd.SetTag(12));
    check_nothrow(bnd.GetMolecule());
    check_nothrow(bnd.GetOrder());
    check_nothrow(bnd.SetOrder(BondOrder::SINGLE));
    check_nothrow(bnd.GetSourceAtom());
    check_nothrow(bnd.GetTargetAtom());
    check_nothrow(bnd.GetAromaticity());
    check_nothrow(bnd.SetAromaticity(true));
    check_nothrow(bnd.GetStereochemistry());
    check_nothrow(bnd.SetStereochemistry(BondStereo::NONE));
    check_nothrow(bnd.SwapSourceTarget());
    check_nothrow(bnd.GetAtoms());
    check_nothrow(bnd.NumAtoms());
    
    // Check correctness of gets
    check_eq(12, bnd.GetTag());
    check_eq(12, bnd.get_tag());
    check_eq(mol, bnd.GetMolecule());
    check_eq(mol, bnd.get_mol().lock());
    check_eq(BondOrder::SINGLE, bnd.GetOrder());
    check_eq(BondOrder::SINGLE, bnd.get_order());
    check_eq(b, bnd.GetSourceAtom());
    check_eq(a, bnd.GetTargetAtom());
    check_eq(a, bnd.get_atms()[1].lock());
    check_eq(b, bnd.get_atms()[0].lock());
    check_eq(true, bnd.GetAromaticity());
    check_eq(true, bnd.get_aromatic());
    check_eq(BondStereo::NONE, bnd.GetStereochemistry());
    check_eq(BondStereo::NONE, bnd.get_stereo());
    check_eq(std::make_pair(b,a), bnd.GetAtoms());
    
    // Check no owning of molecule
    mol.reset();
    check(bnd.get_mol().expired());
    check_eq(Molecule(), bnd.GetMolecule());
    
    // Check num atoms doesn't depend on being active
    check_eq(2, bnd.NumAtoms());
    a.reset(); b.reset();
    check_eq(2, bnd.NumAtoms());
  }
  
  std::string IXBond::ToString() const {
    std::stringstream ss;
    Atom a, b;
    std::tie(a,b) = GetAtoms();
    ss << "Bond(";
    if (!a || !b) ss << "MALFORMED";
    else ss << a->ToString() << ", " << b->ToString();
    ss << ")";
    return ss.str();
  }
  
  std::ostream& operator<<(std::ostream& os, const IXBond& bnd) {
    auto atoms = bnd.GetAtoms();
    os << "Bond(";
    if (atoms.first) os << atoms.first->GetIndex();
    os << ", ";
    if (atoms.second) os << atoms.second->GetIndex();
    os << ")";
    return os;
  }
  
  test_case_fixture(test::BondTestFixture, "IXBond printing methods") {
    // Check ordering is correct
    std::stringstream ss; ss << bnd.imp;
    check_eq("Bond(Atom(0, C), Atom(1, O))", bnd.ToString());
    check_eq("Bond(0, 1)", ss.str());
    bnd.SwapSourceTarget();
    ss.str(""); ss << bnd.imp;
    check_eq("Bond(Atom(1, O), Atom(0, C))", bnd.ToString());
    check_eq("Bond(1, 0)", ss.str());
    
    // Check having bad atoms is handled correctly
    a.reset();
    ss.str(""); ss << bnd.imp;
    check_eq("Bond(MALFORMED)", bnd.ToString());
    check_eq("Bond(1, )", ss.str());
    
    // Check empty bond to ostream does nothing
    ss.str("");
    check_nothrow(ss << Bond());
    check_eq("", ss.str());
  }
  
  void IXBond::Clear() {
    _mol.reset();
    _tag = 0;
    _order = Order::UNDEFINED;
    _stereo = Stereo::UNDEFINED;
    _aromatic = false;
    _atms.fill(_Atom());
  }
  
  test_case_fixture(test::BondTestFixture, "IXBond clearing methods") {
    bnd.SetTag(12);
    bnd.SetOrder(BondOrder::DOUBLE);
    bnd.SetAromaticity(true);
    bnd.SetStereochemistry(BondStereo::E);
    // Pre checks
    check_eq(mol, bnd.get_mol().lock());
    check_eq(12, bnd.get_tag());
    check_eq(BondOrder::DOUBLE, bnd.get_order());
    check_eq(true, bnd.get_aromatic());
    check_eq(BondStereo::E, bnd.get_stereo());
    check_eq(a, bnd.get_atms()[0].lock());
    check_eq(b, bnd.get_atms()[1].lock());
    
    bnd.Clear();
    // Post checks
    check_ne(mol, bnd.get_mol().lock());
    check_ne(12, bnd.get_tag());
    check_ne(BondOrder::DOUBLE, bnd.get_order());
    check_ne(true, bnd.get_aromatic());
    check_ne(BondStereo::E, bnd.get_stereo());
    check_ne(a, bnd.get_atms()[0].lock());
    check_ne(b, bnd.get_atms()[1].lock());
  }
  
  test_suite_close();
  
}

