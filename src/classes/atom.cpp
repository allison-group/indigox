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
#include <indigox/utils/serialise.hpp>

namespace indigox {
  
  template<typename Archive>
  void IXAtom::Serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("molecule", _mol),
            INDIGOX_SERIAL_NVP("formal_charge", _fc),
            INDIGOX_SERIAL_NVP("tag", _tag),
            INDIGOX_SERIAL_NVP("implicit_h_count", _implicitH),
            INDIGOX_SERIAL_NVP("name", _name),
            INDIGOX_SERIAL_NVP("position_x", _pos.x),
            INDIGOX_SERIAL_NVP("position_y", _pos.y),
            INDIGOX_SERIAL_NVP("position_z", _pos.z),
            INDIGOX_SERIAL_NVP("partial_charge", _partial),
            INDIGOX_SERIAL_NVP("stereochemistry", _stereo),
            INDIGOX_SERIAL_NVP("is_aromatic", _aromatic)
            );
    
    // Things that are different between loading and saving
    string_ element;
    if (INDIGOX_IS_OUTPUT_ARCHIVE) element = GetElement()->GetName();
    archive(INDIGOX_SERIAL_NVP("element", element));
    if (INDIGOX_IS_INPUT_ARCHIVE) SetElement(element);
  
    std::vector<uint_> tags; tags.reserve(_dhds.size());
    Molecule m = _mol.lock();
    
    // Bonds
    if (INDIGOX_IS_OUTPUT_ARCHIVE) {
      for (auto it = GetBondIters(); it.first != it.second; ++it.first)
        tags.emplace_back(it.first->lock()->GetTag());
    }
    archive(INDIGOX_SERIAL_NVP("bonds", tags));
    if (INDIGOX_IS_INPUT_ARCHIVE) {
      _bnds.reserve(tags.size());
      for (uint_ t : tags) _bnds.emplace_back(m->GetBond(t));
    }
    tags.clear();
    
    // Angles
    if (INDIGOX_IS_OUTPUT_ARCHIVE) {
      for (auto it = GetAngleIters(); it.first != it.second; ++it.first)
        tags.emplace_back(it.first->lock()->GetTag());
    }
    archive(INDIGOX_SERIAL_NVP("angles", tags));
    if (INDIGOX_IS_INPUT_ARCHIVE) {
      _angs.reserve(tags.size());
      for (uint_ t : tags) _angs.emplace_back(m->GetAngle(t));
    }
    tags.clear();
    
    // Dihedrals
    if (INDIGOX_IS_OUTPUT_ARCHIVE) {
      for (auto it = GetDihedralIters(); it.first != it.second; ++it.first)
        tags.emplace_back(it.first->lock()->GetTag());
    }
    archive(INDIGOX_SERIAL_NVP("dihedrals", tags));
    if (INDIGOX_IS_INPUT_ARCHIVE) {
      _dhds.reserve(tags.size());
      for (uint_ t : tags) _dhds.emplace_back(m->GetDihedral(t));
    }
    tags.clear();
    
  }
  
  INDIGOX_SERIALISE(IXAtom);
  
  inline void __set_property_modified(_Molecule mol, MolProperty p) {
    if (!mol.expired()) mol.lock()->SetPropertyModified(p);
  }
  
  IXAtom::IXAtom(Molecule m) : utils::IXCountableObject<IXAtom>(), _mol(m),
  _elem(), _fc(0), _tag(0), _implicitH(0), _name(), _pos({0.0,0.0,0.0}),
  _partial(0.0), _stereo(Stereo::UNDEFINED), _aromatic(false) { }
  
  string_ IXAtom::ToString() {
    std::stringstream ss;
    ss << "Atom(" << GetIndex() << ", " << GetElement()->GetSymbol() << ")";
    return ss.str();
  }

  void IXAtom::SetElement(Element e) {
    if (e != GetElement()) {
      _elem = e;
      __set_property_modified(_mol, MolProperty::ATOM_ELEMENTS);
    }
  }
  
  size_ IXAtom::GetIndex() const {
    Molecule mol = _mol.lock();
    if (!mol) return GetTag();
    auto be = mol->GetAtoms();
    auto pos = std::find(be.first, be.second, shared_from_this());
    if (pos == be.second) return GetTag();
    return std::distance(be.first, pos);
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
    if (atom) os << "Atom(" << atom->GetIndex() << ")";
    return os;
  }
}
