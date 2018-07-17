#include <sstream>
#include <stdexcept>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/serialise.hpp>

namespace indigox {
  using namespace indigox::utils; // for IXCountableObject and WeakContains
  
  template <typename Archive>
  void IXBond::Serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("molecule", _mol),
            INDIGOX_SERIAL_NVP("tag", _tag),
            INDIGOX_SERIAL_NVP("bond_order", _order),
            INDIGOX_SERIAL_NVP("is_aromatic", _aromatic),
            INDIGOX_SERIAL_NVP("stereochemistry", _stereo)
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
    }
  }
  
  INDIGOX_SERIALISE(IXBond);
  
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

