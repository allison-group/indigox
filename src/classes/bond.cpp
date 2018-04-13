#include <sstream>

#include "indigox/classes/atom.hpp"
#include "indigox/classes/bond.hpp"
#include "indigox/utils/counter.hpp"

namespace indigox {
  
  /*
   *  Bond implementation
   */
  
  // Initalisation methods
  IXBond::IXBond() : utils::CountableObject<IXBond>(), _mol(), _idx(0),
                    _order(Order::UNDEFINED), _aromatic(false),
                    _stereo(Stereo::UNDEFINED), _atoms({{_Atom(),_Atom()}})  { }
  
  IXBond::IXBond(Atom a, Atom b) : IXBond() { _atoms = {{a,b}}; }
  
  std::string IXBond::ToString() const {
    std::stringstream ss;
    ss << "Bond(";
    if (_atoms[0].expired() || _atoms[1].expired()) ss << "MALFORMED";
    else
      ss << _atoms[0].lock()->GetName() << ", " << _atoms[1].lock()->GetName();
    ss << ")";
    return ss.str();
  }
  
  void IXBond::Clear() {
    _mol.reset();
    _idx = 0;
    _order = Order::UNDEFINED;
    _aromatic = false;
    _atoms.fill(_Atom());
    _angles.clear();
    _dihedrals.clear();
  }
  
  void IXBond::RemoveAngle(Angle a) {
    BondAngleIter it = _angles.begin();
    for (; it != _angles.end(); ++it) {
      if (!it->expired() && it->lock() == a) break;
    }
    if (it != _angles.end()) _angles.erase(it);
  }
  
  void IXBond::RemoveDihedral(Dihedral d) {
    BondDihedralIter it = _dihedrals.begin();
    for (; it != _dihedrals.end(); ++it) {
      if (!it->expired() && it->lock() == d) break;
    }
    if (it != _dihedrals.end()) _dihedrals.erase(it);
  }
  
}

