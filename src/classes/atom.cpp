#include <algorithm>
#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>

using namespace indigox;

IXAtom::IXAtom() : utils::CountableObject<IXAtom>(), _mol(), _elem(), _fc(0),
                   _tag(0), _implicitH(0), _name("ATOM"), _pos({0.0,0.0,0.0}),
                   _partial(0.0), _stereo(Stereo::UNDEFINED), _aromatic(false)
                  { }

IXAtom::IXAtom(Molecule m) : IXAtom() {
  _mol = m;
}

string_ IXAtom::ToString() {
  std::stringstream ss;
  Element e = GetElement();
  ss << "Atom(" << _name << ", " << e->GetSymbol() << ")";
  return ss.str();
}


void IXAtom::RemoveBond(Bond b) {
  AtomBondIter it = _bonds.begin();
  for (; it != _bonds.end(); ++it) {
    if (!it->expired() && it->lock() == b) break;
  }
  if (it != _bonds.end()) _bonds.erase(it);
}

void IXAtom::RemoveAngle(Angle a) {
  AtomAngleIter it = _angles.begin();
  for (; it != _angles.end(); ++it) {
    if (!it->expired() && it->lock() == a) break;
  }
  if (it != _angles.end()) _angles.erase(it);
}

void IXAtom::RemoveDihedral(Dihedral d) {
  AtomDihedralIter it = _dihedrals.begin();
  for (; it != _dihedrals.end(); ++it) {
    if (!it->expired() && it->lock() == d) break;
  }
  if (it != _dihedrals.end()) _dihedrals.erase(it);
}

void IXAtom::SetElement(uint_ e)  {
  _elem = IXPeriodicTable::GetInstance()->GetElement(e);
}

void IXAtom::SetElement(string_ e) {
  _elem = IXPeriodicTable::GetInstance()->GetElement(e);
}

Element IXAtom::GetElement() const {
  if(_elem.expired())
    return IXPeriodicTable::GetInstance()->GetUndefinedElement();
  return _elem.lock();
}

void IXAtom::Clear() {
  _mol.reset();
  _elem.reset();
  _fc = 0;
  _tag = 0;
  _implicitH = 0;
  _name = "ATOM";
  _pos = Vec3();
  _partial = 0.0;
  _stereo = Stereo::UNDEFINED;
  _aromatic = false;
  _bonds.clear();
  _angles.clear();
  _dihedrals.clear();
}

void IXAtom::Cleanup() {
  for (auto it = _bonds.begin(); it != _bonds.end(); ) {
    if (it->expired()) it = _bonds.erase(it);
    else ++it;
  }
  for (auto it = _angles.begin(); it != _angles.end(); ) {
    if (it->expired()) it = _angles.erase(it);
    else ++it;
  }
  for (auto it = _dihedrals.begin(); it != _dihedrals.end(); ) {
    if (it->expired()) it = _dihedrals.erase(it);
    else ++it;
  }
}


