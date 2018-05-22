#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox {
  
  IXAtom::IXAtom(Molecule m) : utils::IXCountableObject<IXAtom>(), _mol(m),
  _elem(), _fc(0), _tag(0), _implicitH(0), _name(), _pos({0.0,0.0,0.0}),
  _partial(0.0), _stereo(Stereo::UNDEFINED), _aromatic(false) { }
  
  string_ IXAtom::ToString() {
    std::stringstream ss;
    ss << "Atom(" << _name << ", " << GetElement()->GetSymbol() << ")";
    return ss.str();
  }

  
  void IXAtom::Clear() {
    _mol.reset();
    _elem.reset();
    _fc = 0;
    _tag = 0;
    _implicitH = 0;
    _name = "";
    _pos = Vec3();
    _partial = 0.0;
    _stereo = Stereo::UNDEFINED;
    _aromatic = false;
    _bnds.clear();
    _angs.clear();
    _dhds.clear();
  }
  
  std::ostream& operator<<(std::ostream& os, Atom atom) {
    if (atom) {
      os << "Atom(" << atom->GetUniqueID() << ", "
      << atom->GetElement()->GetSymbol() << ")";
    }
    return os;
  }
}
