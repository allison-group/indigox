#include <sstream>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/utils/serialise.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/dihedral_test.hpp>
#include <indigox/test/forcefield_test.hpp>

namespace indigox {
  
//  test::DihedralTestFixture::DihedralTestFixture()
//  : fftype(test::CreateGenericTestFFDihedral().imp)  {
//    a->SetTag(0); b->SetTag(1); c->SetTag(2); d->SetTag(3);
//    a->SetElement("C"); b->SetElement("O");
//    c->SetElement("F"); d->SetElement("N");
//  }
  
  test_suite_open("Dihedral");
  
  template <typename Archive>
  void Dihedral::save(Archive &archive, const uint32_t) const {
    std::vector<sAtom> atoms;
    for (wAtom at : _atms) atoms.emplace_back(at.lock());
    archive(INDIGOX_SERIAL_NVP("molecule", _mol),
            INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("tag", _tag),
            INDIGOX_SERIAL_NVP("type", _types)
            );
  }
  
  template <typename Archive>
  void Dihedral::load_and_construct(Archive &archive,
                                      cereal::construct<Dihedral> &construct,
                                      const uint32_t) {
    sMolecule m;
    std::vector<sAtom> atoms;
    archive(INDIGOX_SERIAL_NVP("molecule", m),
            INDIGOX_SERIAL_NVP("atoms", atoms));
    
    construct(*atoms[0], *atoms[1], *atoms[2], *atoms[3], *m);
    
    archive(INDIGOX_SERIAL_NVP("tag", construct->_tag),
            INDIGOX_SERIAL_NVP("type", construct->_types));
  }
  INDIGOX_SERIALISE_CONSTRUCT(Dihedral);
  
/*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXDihedral serialisation", T, ixdihed_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::DihedralTestFixture fixture;
    Dihedral saved = fixture.dhd.imp;

    saved->SetTag(45);
    saved->SetType(fixture.fftype);
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved, saved->GetAtoms(), saved->GetMolecule(),
                        fixture.fftype));
    }

    Dihedral loaded;
    test::TestDihedral::Atoms atoms;
    Molecule mol;
    FFDihedral ff_loaded;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, atoms, mol, ff_loaded));
    }
    fixture.dhd.imp = loaded;
    check(!fixture.dhd.get_atms()[0].expired());
    check(!fixture.dhd.get_atms()[1].expired());
    check(!fixture.dhd.get_atms()[2].expired());
    check(!fixture.dhd.get_atms()[3].expired());
    check(!fixture.dhd.get_mol().expired());
    check_eq(saved->GetTag(), loaded->GetTag());
    check_eq(atoms.first->GetTag(), saved->GetAtoms().first->GetTag());
    check_eq(atoms.second->GetTag(), saved->GetAtoms().second->GetTag());
    check_eq(atoms.third->GetTag(), saved->GetAtoms().third->GetTag());
    check_eq(atoms.fourth->GetTag(), saved->GetAtoms().fourth->GetTag());
    check_eq(atoms, loaded->GetAtoms());
    check_eq(ff_loaded, loaded->GetType());
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixdihed_serial, ixserial<IXDihedral>);
 */
  
  Dihedral::Dihedral(Atom& a, Atom& b, Atom& c, Atom& d, Molecule& m) :
  _mol(m.weak_from_this()), _tag(0),
  _atms({a.weak_from_this(), b.weak_from_this(), c.weak_from_this(),
    d.weak_from_this()}) { }
  
/*  test_case_fixture(test::DihedralTestFixture, "IXDihedral construction") {
    check_nothrow(test::TestDihedral tdhd(a,b,c,d,mol));
    test::TestDihedral tdhd(a,b,c,d,mol);
    check_eq(mol, tdhd.get_mol().lock());
    check_eq(0, tdhd.get_tag());
    check_eq(a, tdhd.get_atms()[0].lock());
    check_eq(b, tdhd.get_atms()[1].lock());
    check_eq(c, tdhd.get_atms()[2].lock());
    check_eq(d, tdhd.get_atms()[3].lock());
    check_eq(FFDihedral(), tdhd.get_type());

    // Check unique IDs correctly update
    test::TestDihedral dhd1(a,b,c,d,mol);
    test::TestDihedral dhd2(a,b,c,d,mol);
    check_ne(dhd1.GetUniqueID(), dhd2.GetUniqueID());
    check_eq(dhd1.GetUniqueID() + 1, dhd2.GetUniqueID());
  }
 */
  
  size_t Dihedral::GetIndex() const {
    if (_mol.expired()) return GetTag();
    Molecule& mol = *_mol.lock();
    auto be = mol.GetDihedralIters();
    auto pos = std::find(be.first, be.second, shared_from_this());
    if (pos == be.second) return GetTag();
    return std::distance(be.first, pos);
  }
  
/*  test_case_fixture(test::DihedralTestFixture, "IXDihedral getting and setting") {
    // check no throwing
    check_nothrow(dhd.GetTag());
    check_nothrow(dhd.SetTag(72));
    check_nothrow(dhd.GetMolecule());
    check_nothrow(dhd.GetAtoms());
    check_nothrow(dhd.SwapOrder());
    check_nothrow(dhd.NumAtoms());
    check_nothrow(dhd.SetType(fftype));

    // Check correctness of gets
    check_eq(72, dhd.GetTag());
    check_eq(72, dhd.get_tag());
    check_eq(mol, dhd.GetMolecule());
    check_eq(mol, dhd.get_mol().lock());
    check_eq(stdx::make_quad(d, c, b, a), dhd.GetAtoms());
    check_eq(dhd.GetTag(), dhd.GetIndex());
    check_eq(fftype, dhd.GetType());

    // Check correct no owning of molecule
    mol.reset();
    check(dhd.get_mol().expired());
    check_eq(Molecule(), dhd.GetMolecule());

    // Check still gets tag when mol is dead
    check_eq(dhd.GetTag(), dhd.GetIndex());

    // Check num atoms doesn't depend on being active
    check_eq(4, dhd.NumAtoms());
    a.reset(); b.reset(); c.reset(); d.reset();
    check_eq(4, dhd.NumAtoms());
  }
 */
  
  std::string Dihedral::ToString() const {
    std::stringstream ss;
    Atom& a = *_atms[0].lock();
    Atom& b = *_atms[1].lock();
    Atom& c = *_atms[2].lock();
    Atom& d = *_atms[3].lock();
    ss << "Dihedral(" << a.ToString() << ", " << b.ToString() << ", "
    << c.ToString() << ", " << d.ToString() << ")";
    return ss.str();
  }
  
  std::ostream& operator<<(std::ostream& os, const Dihedral& dhd) {
    auto s = dhd.GetAtoms();
    os << "Dihedral(" << s.first.GetIndex() << ", " << s.second.GetIndex()
       << ", " << s.third.GetIndex() << ", " << s.fourth.GetIndex() << ")";
    return os;
  }
  
/*  test_case_fixture(test::DihedralTestFixture, "IXDihedral printing methods") {
    // Check ordering is correct
    std::stringstream ss; ss << dhd.imp;
    check_eq("Dihedral(Atom(0, C), Atom(1, O), Atom(2, F), Atom(3, N))",
             dhd.ToString());
    check_eq("Dihedral(0, 1, 2, 3)", ss.str());
    dhd.SwapOrder();
    ss.str(""); ss << dhd.imp;
    check_eq("Dihedral(Atom(3, N), Atom(2, F), Atom(1, O), Atom(0, C))",
             dhd.ToString());
    check_eq("Dihedral(3, 2, 1, 0)", ss.str());

    // Check having bad atoms is handled correctly
    a.reset(); d.reset();
    ss.str(""); ss << dhd.imp;
    check_eq("Dihedral(MALFORMED)", dhd.ToString());
    check_eq("Dihedral(, 2, 1, )", ss.str());

    // Check empty dihedral to ostream does nothing
    ss.str("");
    check_nothrow(ss << Dihedral());
    check_eq("", ss.str());
  }
 */
  
  FFDihedral& Dihedral::GetType(size_t pos) const {
    if (pos >= _types.size()) throw std::out_of_range("Invalid position");
    return *_types.at(pos).lock();
  }
  
  void Dihedral::AddType(FFDihedral& type) {
    if (!GetMolecule().HasForcefield())
      GetMolecule().SetForcefield(type.GetForcefield());
    if (type.GetForcefield().shared_from_this()
        != GetMolecule().GetForcefield().shared_from_this())
      throw std::runtime_error("Dihedral type is not from molecule's forcefield");
    
    _types.emplace_back(type.weak_from_this());
    // Keep the list of dihedrals ordered
    std::sort(_types.begin(), _types.end(),
              [](wFFDihedral a, wFFDihedral b){
                FFDihedral& ar = *a.lock();
                FFDihedral& br = *b.lock();
                if (ar.GetType() != br.GetType())
                  return ar.GetType() < br.GetType();
                return ar.GetID() < br.GetID();
              });
  }
  
  void Dihedral::RemoveType(FFDihedral& type) {
    auto pos = utils::WeakContainsShared(_types.begin(), _types.end(),
                                         type.shared_from_this());
    if (pos != _types.end()) _types.erase(pos);
  }
  
  void Dihedral::Clear() {
    _mol.reset();
    _atms.fill(wAtom());
    _types.clear();
  }
  
  test_suite_close();
}
