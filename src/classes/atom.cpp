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
  using namespace utils;  // for WeakContainsShared(...)
  
//  std::ostream& operator<<(std::ostream& os, IXAtom::Stereo s) {
//    switch (s) {
//      case IXAtom::Stereo::UNDEFINED:
//        return (os << "UNDEFINED");
//      case IXAtom::Stereo::ACHIRAL:
//        return (os << "ACHIRAL");
//      case IXAtom::Stereo::R:
//        return (os << "R");
//      case IXAtom::Stereo::S:
//        return (os << "S");
//    }
//  }
  
  IXAtom::IXAtom() : utils::IXCountableObject<IXAtom>(), _mol(), _elem(),
  _fc(0), _tag(0), _implicitH(0), _name(), _pos({0.0,0.0,0.0}), _partial(0.0),
  _stereo(Stereo::UNDEFINED), _aromatic(false) { }
  
  IXAtom::IXAtom(Molecule m) : utils::IXCountableObject<IXAtom>(), _mol(m),
  _elem(), _fc(0), _tag(0), _implicitH(0), _name(), _pos({0.0,0.0,0.0}),
  _partial(0.0), _stereo(Stereo::UNDEFINED), _aromatic(false) { }
  
  string_ IXAtom::ToString() {
    std::stringstream ss;
    ss << "Atom(" << _name << ", " << GetElement()->GetSymbol() << ")";
    return ss.str();
  }
  
  bool IXAtom::AddBond(Bond b) {
    if (!b) return false;
    AtomBondIter it = WeakContainsShared(_bnds.begin(), _bnds.end(), b);
    if (it == _bnds.end()) return static_cast<void>(_bnds.push_back(b)), true;
    return false;
  }
  
  bool IXAtom::RemoveBond(Bond b) {
    if (!b) return false;
    AtomBondIter it = WeakContainsShared(_bnds.begin(), _bnds.end(), b);
    if (it != _bnds.end()) return static_cast<void>(_bnds.erase(it)), true;
    return false;
  }
  
  size_ IXAtom::NumBonds() const {
    return std::count_if(_bnds.begin(), _bnds.end(),
                         [](_Bond b) { return !b.expired(); });
  }
  
  bool IXAtom::AddAngle(Angle a) {
    if (!a) return false;
    AtomAngleIter it = WeakContainsShared(_angs.begin(), _angs.end(), a);
    if (it == _angs.end()) return static_cast<void>(_angs.push_back(a)), true;
    return false;
  }
  
  bool IXAtom::RemoveAngle(Angle a) {
    if (!a) return false;
    AtomAngleIter it = WeakContainsShared(_angs.begin(), _angs.end(), a);
    if (it != _angs.end()) return static_cast<void>(_angs.erase(it)), true;
    return false;
  }
  
  size_ IXAtom::NumAngles() const {
    return std::count_if(_angs.begin(), _angs.end(),
                         [](_Angle a) { return !a.expired(); });
  }
  
  bool IXAtom::AddDihedral(Dihedral d) {
    if (!d) return false;
    AtomDihedralIter it = WeakContainsShared(_dhds.begin(), _dhds.end(), d);
    if (it == _dhds.end()) return static_cast<void>(_dhds.push_back(d)), true;
    return false;
  }
  
  bool IXAtom::RemoveDihedral(Dihedral d) {
    AtomDihedralIter it = WeakContainsShared(_dhds.begin(), _dhds.end(), d);
    if (it != _dhds.end()) return static_cast<void>(_dhds.erase(it)), true;
    return false;
  }
  
  size_ IXAtom::NumDihedrals() const {
    return std::count_if(_dhds.begin(), _dhds.end(),
                         [](_Dihedral d) { return !d.expired(); });
  }
  
  void IXAtom::SetElement(Element e) {
    _elem = e;
  }
  
  void IXAtom::SetElement(uint_ e)  {
    _elem = GetPeriodicTable()->GetElement(e);
  }
  
  void IXAtom::SetElement(string_ e) {
    _elem = GetPeriodicTable()->GetElement(e);
  }
  
  Element IXAtom::GetElement() const {
    if(_elem.expired())
      return GetPeriodicTable()->GetUndefinedElement();
    return _elem.lock();
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
  
  void IXAtom::Cleanup() {
    _bnds.erase(std::remove_if(_bnds.begin(), _bnds.end(),
                                [](_Bond b){ return b.expired(); }),
                 _bnds.end());
    _angs.erase(std::remove_if(_angs.begin(), _angs.end(),
                                 [](_Angle a){ return a.expired(); }),
                  _angs.end());
    _dhds.erase(std::remove_if(_dhds.begin(), _dhds.end(),
                                    [](_Dihedral d){ return d.expired(); }),
                     _dhds.end());
  }
  
  std::ostream& operator<<(std::ostream& os, Atom atom) {
    if (atom) {
      os << "Atom(" << atom->GetUniqueID() << ", "
      << atom->GetElement()->GetSymbol() << ")";
    }
    return os;
  }
}
