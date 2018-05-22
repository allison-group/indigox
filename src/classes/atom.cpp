#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox {
  
  inline void __set_property_modified(_Molecule mol, IXMolecule::Property p) {
    if (!mol.expired()) mol.lock()->SetPropertyModified(p);
  }
  
  IXAtom::IXAtom(Molecule m) : utils::IXCountableObject<IXAtom>(), _mol(m),
  _elem(), _fc(0), _tag(0), _implicitH(0), _name(), _pos({0.0,0.0,0.0}),
  _partial(0.0), _stereo(Stereo::UNDEFINED), _aromatic(false) { }
  
  string_ IXAtom::ToString() {
    std::stringstream ss;
    ss << "Atom(" << _name << ", " << GetElement()->GetSymbol() << ")";
    return ss.str();
  }

  void IXAtom::SetElement(Element e) {
    if (e != GetElement()) {
      _elem = e;
      __set_property_modified(_mol, IXMolecule::Property::ATOM_ELEMENTS);
    }
  }
  
  void IXAtom::SetElement(string_ e) {
    Element elem = GetPeriodicTable()->GetElement(e);
    if (e != GetElement()) {
      _elem = elem;
      __set_property_modified(_mol, IXMolecule::Property::ATOM_ELEMENTS);
    }
  }
  
  void IXAtom::SetElement(uint_ e) {
    Element elem = GetPeriodicTable()->GetElement(e);
    if (e != GetElement()) {
      _elem = elem;
      __set_property_modified(_mol, IXMolecule::Property::ATOM_ELEMENTS);
    }
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
