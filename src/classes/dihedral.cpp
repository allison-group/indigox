#include <sstream>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox {
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
