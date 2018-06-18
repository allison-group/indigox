#include <sstream>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox {
  IXAngle::IXAngle(Atom a, Atom b, Atom c, Molecule m)
  : utils::IXCountableObject<IXAngle>(), _mol(m), _tag(0), _atms({{a,b,c}}) { }
  
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
  
  void IXAngle::Clear() {
    _mol.reset();
    _tag = 0;
    _atms[0].reset();
    _atms[1].reset();
    _atms[2].reset();
  }
}
