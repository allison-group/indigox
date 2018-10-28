#include <sstream>
#include <stdexcept>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/bond_test.hpp>
#include <indigox/test/forcefield_test.hpp>

namespace indigox {
//  test::BondTestFixture::BondTestFixture()
//  : fftype(test::CreateGenericTestFFBond().imp) {
//    a->SetTag(0); b->SetTag(1);
//    a->SetElement("C"); b->SetElement("O");
//  }
  
  using namespace indigox::utils; // for IXCountableObject and WeakContains
  
  test_suite_open("Bond");
  
  // Serialisation of IXBond
  template <typename Archive>
  void Bond::save(Archive &archive, const uint32_t) const {
    std::vector<sAtom> atoms;
    for (wAtom at : _atms) atoms.emplace_back(at.lock());
    
    archive(INDIGOX_SERIAL_NVP("molecule", _mol.lock()),
            INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("tag", _tag),
            INDIGOX_SERIAL_NVP("bond_order", _order),
            INDIGOX_SERIAL_NVP("is_aromatic", _aromatic),
            INDIGOX_SERIAL_NVP("stereochemistry", _stereo),
            INDIGOX_SERIAL_NVP("type", _type)
            );
  }
  
  template <typename Archive>
  void Bond::load_and_construct(Archive &archive,
                                  cereal::construct<Bond> &construct,
                                  const uint32_t) {
    sMolecule m;
    std::vector<sAtom> atoms;
    archive(INDIGOX_SERIAL_NVP("molecule", m),
            INDIGOX_SERIAL_NVP("atoms", atoms));
    
    construct(*atoms[0], *atoms[1], *m);
    archive(INDIGOX_SERIAL_NVP("tag", construct->_tag),
            INDIGOX_SERIAL_NVP("bond_order", construct->_order),
            INDIGOX_SERIAL_NVP("is_aromatic", construct->_aromatic),
            INDIGOX_SERIAL_NVP("stereochemistry", construct->_stereo),
            INDIGOX_SERIAL_NVP("type", construct->_type)
            );
  }
  INDIGOX_SERIALISE_CONSTRUCT(Bond);
  
/*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXBond serialisation", T, ixbond_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::BondTestFixture fixture;
    Bond saved = fixture.bnd.imp;

    saved->SetTag(12);
    saved->SetOrder(BondOrder::TRIPLE);
    saved->SetAromaticity(true);
    saved->SetStereochemistry(BondStereo::Z);
    saved->SetType(fixture.fftype);

    std::ostringstream os;
    {
      Out oar(os);
      // atoms and molecule saved so they can be kept alive on load
      check_nothrow(oar(saved, saved->GetSourceAtom(), saved->GetTargetAtom(),
                        saved->GetMolecule(), fixture.fftype));
    }

    Bond loaded;
    Atom source, target;
    Molecule mol;
    FFBond typ_loaded;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, source, target, mol, typ_loaded));
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
    check_eq(source, loaded->GetSourceAtom());
    check_eq(target, loaded->GetTargetAtom());
    check_eq(mol, loaded->GetMolecule());
    check_eq(typ_loaded, loaded->GetType());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixbond_serial, ixserial<IXBond>);
 */
  
  Bond::Bond(Atom& a, Atom& b, Molecule& m) :
  IXCountableObject<Bond>(), _mol(m.weak_from_this()), _tag(0),
  _order(Order::UNDEFINED), _aromatic(false), _stereo(Stereo::UNDEFINED),
  _atms({a.weak_from_this(),b.weak_from_this()}) { }
  
/*  test_case_fixture(test::BondTestFixture, "IXBond construction") {
    check_nothrow(test::TestBond tbnd(a,b,mol));
    test::TestBond tbnd(a,b,mol);
    check_eq(mol, tbnd.get_mol().lock());
    check_eq(0, tbnd.get_tag());
    check_eq(BondOrder::UNDEFINED, tbnd.get_order());
    check_eq(false, tbnd.get_aromatic());
    check_eq(BondStereo::UNDEFINED, tbnd.get_stereo());
    check_eq(a, tbnd.get_atms()[0].lock());
    check_eq(b, tbnd.get_atms()[1].lock());
    check_eq(FFBond(), tbnd.get_type());

    // Check unique IDs correctly update
    test::TestBond bnd1(a,b,mol);
    test::TestBond bnd2(a,b,mol);
    check_ne(bnd1.GetUniqueID(), bnd2.GetUniqueID());
    check_eq(bnd1.GetUniqueID() + 1, bnd2.GetUniqueID());
  }
 */
  
  size_t Bond::GetIndex() const {
    if (!_mol.expired()) return GetTag();
    Molecule& mol = *_mol.lock();
    auto be = mol.GetBondIters();
    auto pos = std::find(be.first, be.second, shared_from_this());
    if (pos == be.second) return GetTag();
    return std::distance(be.first, pos);
  }
  
/*  test_case_fixture(test::BondTestFixture, "IXBond getting and setting") {
    // Check no throwing
    check_nothrow(bnd.SetTag(12));
    check_nothrow(bnd.SetOrder(BondOrder::SINGLE));
    check_nothrow(bnd.SetAromaticity(true));
    check_nothrow(bnd.SetStereochemistry(BondStereo::NONE));
    check_nothrow(bnd.SwapSourceTarget());
    check_nothrow(bnd.SetType(fftype));

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
    check_eq(bnd.GetTag(), bnd.GetIndex());
    check_eq(fftype, bnd.GetType());

    // Check no owning of molecule
    mol.reset();
    check(bnd.get_mol().expired());
    check_eq(Molecule(), bnd.GetMolecule());
    // Check still gets tag when mol is dead
    check_eq(bnd.GetTag(), bnd.GetIndex());

    // Check num atoms doesn't depend on being active
    check_eq(2, bnd.NumAtoms());
    a.reset(); b.reset();
    check_eq(2, bnd.NumAtoms());
  }
 */
  
  std::string Bond::ToString() const {
    std::stringstream ss;
    Atom& a = *_atms[0].lock();
    Atom& b = *_atms[1].lock();
    ss << "Bond(" << a.ToString() << ", " << b.ToString() << ")";
    return ss.str();
  }
  
  std::ostream& operator<<(std::ostream& os, const Bond& bnd) {
    auto s = bnd.GetAtoms();
    os << "Bond(" << s.first.GetIndex() << ", " << s.second.GetIndex() << ")";
    return os;
  }
  
/*  test_case_fixture(test::BondTestFixture, "IXBond printing methods") {
    // Check ordering is correct
    std::stringstream ss;
    check_nothrow(ss << bnd.imp);
    check_eq("Bond(Atom(0, C), Atom(1, O))", bnd.ToString());
    check_eq("Bond(0, 1)", ss.str());
    bnd.SwapSourceTarget();
    ss.str("");
    check_nothrow(ss << bnd.imp);
    check_eq("Bond(Atom(1, O), Atom(0, C))", bnd.ToString());
    check_eq("Bond(1, 0)", ss.str());

    // Check having bad atoms is handled correctly
    a.reset();
    ss.str("");
    check_nothrow(ss << bnd.imp);
    check_eq("Bond(MALFORMED)", bnd.ToString());
    check_eq("Bond(1, )", ss.str());

    // Check empty bond to ostream does nothing
    ss.str("");
    check_nothrow(ss << Bond());
    check_eq("", ss.str());
  }
 */
  
  void Bond::SetType(FFBond& type) {
    if (!GetMolecule().HasForcefield())
      GetMolecule().SetForcefield(type.GetForcefield());
    if (type.GetForcefield().shared_from_this()
        != GetMolecule().GetForcefield().shared_from_this())
      throw std::runtime_error("Bond type is not from molecule's forcefield");
    _type = type.weak_from_this();
  }
  
  void Bond::Clear() {
    _mol.reset();
    _order = Order::UNDEFINED;
    _stereo = Stereo::UNDEFINED;
    _type.reset();
    _atms.fill(wAtom());
  }
  
  test_suite_close();
  
}

