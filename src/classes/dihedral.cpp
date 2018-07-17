#include <sstream>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/serialise.hpp>

namespace indigox {
  
  template <typename Archive>
  void IXDihedral::Serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("molecule", _mol),
            INDIGOX_SERIAL_NVP("tag", _tag)
            );
    
    // Things that are different between save and load
    std::vector<uint_> atoms;
    if (INDIGOX_IS_OUTPUT_ARCHIVE)
      for (_Atom at : _atms) atoms.emplace_back(at.lock()->GetTag());
    archive(INDIGOX_SERIAL_NVP("atoms", atoms));
    if (INDIGOX_IS_INPUT_ARCHIVE) {
      Molecule m = _mol.lock();
      _atms[0] = m->GetAtom(atoms[0]);
      _atms[1] = m->GetAtom(atoms[1]);
      _atms[2] = m->GetAtom(atoms[2]);
      _atms[3] = m->GetAtom(atoms[3]);
    }
  }
  
  INDIGOX_SERIALISE(IXDihedral);
  
  IXDihedral::IXDihedral(Atom a, Atom b, Atom c, Atom d, Molecule m)
  : _mol(m), _tag(0), _atms({{a,b,c,d}}) { }
  
  string_ IXDihedral::ToString() const {
    std::stringstream ss;
    Atom a = _atms[0].lock();
    Atom b = _atms[1].lock();
    Atom c = _atms[2].lock();
    Atom d = _atms[3].lock();
    ss << "Dihedral(";
    if (!a || !b || !c || !d) ss << "MALFORMED";
    else {
      ss << a->ToString() << ", " << b->ToString() << ", " << c->ToString()
      << ", " << d->ToString();
    }
    ss << ")";
    return ss.str();
  }
  
  void IXDihedral::Clear() {
    _mol.reset();
    _tag = 0;
    _atms[0].reset();
    _atms[1].reset();
    _atms[2].reset();
    _atms[3].reset();
  }
}
