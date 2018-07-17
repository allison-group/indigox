#include <sstream>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/serialise.hpp>

namespace indigox {
  
  template <typename Archive>
  void IXAngle::Serialise(Archive &archive, const uint32_t) {
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
    }
  }
  
  INDIGOX_SERIALISE(IXAngle);
  
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
