#include <sstream>
#include <stdexcept>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/common.hpp>

namespace indigox {
  using namespace indigox::utils; // for IXCountableObject and WeakContains
  
  IXBond::IXBond(Atom a, Atom b, Molecule m) : IXCountableObject<IXBond>(),
  _mol(m), _tag(0), _order(Order::UNDEFINED), _aromatic(false),
  _stereo(Stereo::UNDEFINED), _atms({{a,b}}) { }
  
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
  }
  
  std::ostream& operator<<(std::ostream& os, Bond bnd) {
    if (bnd) {
      os << "Bond(" << bnd->GetSourceAtom() << ", " << bnd->GetTargetAtom();
      os << ")";
    }
    return os;
  }
}

