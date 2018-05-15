#include <sstream>
#include <stdexcept>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/common.hpp>

namespace indigox {
  using namespace indigox::utils; // for IXCountableObject and WeakContains
  
  IXBond::IXBond() : IXCountableObject<IXBond>(), _mol(), _tag(0),
  _order(Order::UNDEFINED), _aromatic(false), _stereo(Stereo::UNDEFINED),
  _atms({{_Atom(),_Atom()}})  { }
  
  IXBond::IXBond(Atom a, Atom b) : IXCountableObject<IXBond>(), _mol(),
  _tag(0), _order(Order::UNDEFINED), _aromatic(false),
  _stereo(Stereo::UNDEFINED), _atms({{a,b}})  {
    if (a && a == b)
      throw std::logic_error("Cannot create a bond with the same atom.");
  }
  
  bool IXBond::SetSourceAtom(Atom atom) {
    if (atom && atom != _atms[1].lock()) {
      _atms[0] = atom;
      return true;
    }
    return false;
  }
  
  bool IXBond::SetTargetAtom(Atom atom) {
    if (atom && atom != _atms[0].lock()) {
      _atms[1] = atom;
      return true;
    }
    return false;
  }
  
  std::string IXBond::ToString() const {
    std::stringstream ss;
    ss << "Bond(";
    if (_atms[0].expired() || _atms[1].expired()) ss << "MALFORMED";
    else
      ss << _atms[0].lock()->ToString() << ", " << _atms[1].lock()->ToString();
    ss << ")";
    return ss.str();
  }
  
  void IXBond::Clear() {
    _mol.reset();
    _tag = 0;
    _order = Order::UNDEFINED;
    _stereo = Stereo::UNDEFINED;
    _aromatic = false;
    _atms.fill(_Atom());
    _angs.clear();
    _dhds.clear();
  }
  
  size_ IXBond::NumAtoms() const {
    return std::count_if(_atms.begin(), _atms.end(),
                         [](_Atom a) { return !a.expired(); });
  }
  
  bool IXBond::AddAngle(Angle a) {
    if (!a) return false;
    BondAngleIter it = WeakContainsShared(_angs.begin(), _angs.end(), a);
    if (it == _angs.end()) return static_cast<void>(_angs.push_back(a)), true;
    return false;
  }
  
  bool IXBond::RemoveAngle(Angle a) {
    if (!a) return false;
    BondAngleIter it = WeakContainsShared(_angs.begin(), _angs.end(), a);
    if (it != _angs.end()) return static_cast<void>(_angs.erase(it)), true;
    return false;
  }
  
  size_ IXBond::NumAngles() const {
    return std::count_if(_angs.begin(), _angs.end(),
                         [](_Angle a) { return !a.expired(); });
  }
  
  bool IXBond::AddDihedral(Dihedral d) {
    if (!d) return false;
    BondDihedralIter it = WeakContainsShared(_dhds.begin(), _dhds.end(), d);
    if (it == _dhds.end()) return static_cast<void>(_dhds.push_back(d)), true;
    return false;
  }
  
  bool IXBond::RemoveDihedral(Dihedral d) {
    if (!d) return false;
    BondDihedralIter it = WeakContainsShared(_dhds.begin(), _dhds.end(), d);
    if (it != _dhds.end()) return static_cast<void>(_dhds.erase(it)), true;
    return false;
  }
  
  size_ IXBond::NumDihedrals() const {
    return std::count_if(_dhds.begin(), _dhds.end(),
                         [](_Dihedral d) { return  !d.expired(); });
  }
  
  void IXBond::Cleanup() {
    _angs.erase(std::remove_if(_angs.begin(), _angs.end(),
                               [](_Angle a){ return a.expired(); }),
                _angs.end());
    _dhds.erase(std::remove_if(_dhds.begin(), _dhds.end(),
                               [](_Dihedral d){ return d.expired(); }),
                _dhds.end());
  }
  
  std::ostream& operator<<(std::ostream& os, Bond bnd) {
    if (bnd) {
      os << "Bond(" << bnd->GetSourceAtom() << ", " << bnd->GetTargetAtom();
      os << ")";
    }
    return os;
  }
}

