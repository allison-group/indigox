#include <sstream>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/serialise.hpp>

#ifndef INDIGOX_DISABLE_TESTS
#include <indigox/test/angle_test.hpp>
#endif

namespace indigox {
  
  // Serialisation of IXAngle
  template <typename Archive>
  void IXAngle::save(Archive &archive, const uint32_t) const {
    // Serialising molecule first sets all the tags correctly
    archive(INDIGOX_SERIAL_NVP("molecule", _mol));
    std::vector<uint_> atoms;
    for (_Atom at : _atms) atoms.emplace_back(at.lock()->GetTag());
    archive(INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("tag", _tag));
  }
  
  template <typename Archive>
  void IXAngle::load(Archive &, const uint32_t) { }
  
  template <typename Archive>
  void IXAngle::load_and_construct(Archive &archive, cereal::construct<IXAngle> &construct, const uint32_t) {
    std::vector<uint_> atoms;
    uid_ t;
    Molecule m;
    
    archive(INDIGOX_SERIAL_NVP("molecule", m),
            INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("tag", t)
            );

    construct(m->GetAtom(atoms[0]),
              m->GetAtom(atoms[1]),
              m->GetAtom(atoms[2]),
              m);
    
    construct->_tag = t;
  }
  
  INDIGOX_SERIALISE_SPLIT(IXAngle);
  
  test_suite_open("IXAngle");
  
  IXAngle::IXAngle(Atom a, Atom b, Atom c, Molecule m)
  : utils::IXCountableObject<IXAngle>(), _mol(m), _tag(0), _atms({{a,b,c}}) { }

  test_case_fixture(test::AngleTestFixture, "IXAngle construction") {
    check_nothrow(test::TestAngle(a,b,c,mol));
    check_eq(a, ang.get_atms()[0].lock());
    check_ne(a, ang.get_atms()[1].lock());
    check_ne(a, ang.get_atms()[2].lock());
    check_eq(b, ang.get_atms()[1].lock());
    check_ne(b, ang.get_atms()[0].lock());
    check_ne(b, ang.get_atms()[2].lock());
    check_eq(c, ang.get_atms()[2].lock());
    check_ne(c, ang.get_atms()[0].lock());
    check_ne(c, ang.get_atms()[1].lock());
    check_eq(0, ang.get_tag());
    check_eq(mol, ang.get_mol().lock());
  }
  
  test_case_fixture(test::AngleTestFixture, "IXAngle getting and setting") {
    check_eq(0, ang.GetTag());
    check_nothrow(ang.SetTag(72));
    check_eq(72, ang.GetTag());
    check_eq(72, ang.get_tag());
    
    check_eq(mol, ang.GetMolecule());
    mol.reset();
    check(!ang.GetMolecule());
    check_eq(Molecule(), ang.GetMolecule());
    
    auto atoms = ang.GetAtoms();
    check_eq(a, atoms.first);
    check_eq(b, atoms.second);
    check_eq(c, atoms.third);
    atoms = ang_2.GetAtoms();
    check_eq(c, atoms.first);
    check_eq(a, atoms.second);
    check_eq(b, atoms.third);
    ang_2.SwapOrder();
    atoms = ang_2.GetAtoms();
    check_eq(b, atoms.first);
    check_eq(a, atoms.second);
    check_eq(c, atoms.third);
    
    check_eq(3, ang.NumAtoms());
  }
  
  string_ IXAngle::ToString() const {
    std::stringstream ss;
    Atom a = _atms[0].lock();
    Atom b = _atms[1].lock();
    Atom c = _atms[2].lock();
    ss << "Angle(";
    if (!a || !b || !c) ss << "MALFORMED";
    else ss << a->ToString() << ", " << b->ToString() << ", " << c->ToString();
    ss << ")";
    return ss.str();
  }
  
  std::ostream& operator<<(std::ostream& os, const IXAngle& ang) {
    auto atoms = ang.GetAtoms();
    os << "Angle(";
    if (atoms.first) os << atoms.first->GetIndex();
    os << ", ";
    if (atoms.second) os << atoms.second->GetIndex();
    os << ", ";
    if (atoms.third) os << atoms.third->GetIndex();
    os << ")";
    return os;
  }
  
  test_case_fixture(test::AngleTestFixture, "IXAngle printing methods") {
    check_eq("Angle(Atom(0, XX), Atom(1, XX), Atom(2, XX))", ang.ToString());
    ang.SwapOrder();
    check_eq("Angle(Atom(2, XX), Atom(1, XX), Atom(0, XX))", ang.ToString());
    
    std::stringstream ss;
    ss << ang_2.imp;
    check_eq("Angle(2, 0, 1)", ss.str());
    ang_2.SwapOrder(); ss.str(""); ss << ang_2.imp;
    check_eq("Angle(1, 0, 2)", ss.str());
    
    a.reset();
    check_eq("Angle(MALFORMED)", ang.ToString());
    b.reset(); c.reset();
    ss.str(""); ss << ang_2.imp;
    check_eq("Angle(, , )", ss.str());
  }
  
  void IXAngle::Clear() {
    _mol.reset();
    _tag = 0;
    _atms[0].reset();
    _atms[1].reset();
    _atms[2].reset();
  }
  
  test_case_fixture(test::AngleTestFixture, "IXAngle clearing methods") {
    ang.SetTag(72);
    check_eq(ang.get_tag(), 72);
    check_eq(ang.get_mol().lock(), mol);
    check_eq(ang.get_atms()[0].lock(), a);
    check_eq(ang.get_atms()[1].lock(), b);
    check_eq(ang.get_atms()[2].lock(), c);
    
    ang.Clear();
    check_ne(ang.get_tag(), 72);
    check_ne(ang.get_mol().lock(), mol);
    check_ne(ang.get_atms()[0].lock(), a);
    check_ne(ang.get_atms()[1].lock(), b);
    check_ne(ang.get_atms()[2].lock(), c);
    check_eq(ang.get_tag(), 0);
    check_eq(ang.get_mol().lock(), Molecule());
    check_eq(ang.get_atms()[0].lock(), Atom());
    check_eq(ang.get_atms()[1].lock(), Atom());
    check_eq(ang.get_atms()[2].lock(), Atom());
  }
  
  test_suite_close();
}
