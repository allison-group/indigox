#include <algorithm>
#include <cstdint>
#include <memory>
#include <sstream>

#include "indigox/classes/atom.hpp"
#include "indigox/classes/periodictable.hpp"
#include "indigox/utils/counter.hpp"

using namespace indigox;

/*
 *  Atom implementation
 */

// Initalisation methods
IXAtom::IXAtom() : utils::CountableObject<IXAtom>(), _mol(), _elem() {
  _idx = GetUniqueID();
}

IXAtom::IXAtom(Molecule m)
: utils::CountableObject<IXAtom>(), _mol(m), _elem(){
  _idx = GetUniqueID();
}

std::string IXAtom::ToString() {
  std::stringstream ss;
  Element e = GetElement();
#ifndef DEBUG
  ss << "Atom(" << _name << ", " << e->GetSymbol() << ")";
#else
  ss << "Atom(" << _name << "-" << _idx << ", " << e->GetSymbol();
  ss << ", " << _pos.x << ", " << _pos.y << ", " << _pos.z << ")";
#endif
  return ss.str();
}

void IXAtom::RemoveBond(Bond b) {
  AtomBondIter it = _bonds.begin();
  for (; it != _bonds.end(); ++it) {
    if (it->lock() == b) break;
  }
  if (it != _bonds.end()) _bonds.erase(it);
}

void IXAtom::RemoveAngle(Angle a) {
  AtomAngleIter it = _angles.begin();
  for (; it != _angles.end(); ++it) {
    if (it->lock() == a) break;
  }
  if (it != _angles.end()) _angles.erase(it);
}

void IXAtom::Clear() {
  _mol.reset();
  _elem = IXPeriodicTable::GetInstance()->GetElement(0);
  _bonds.clear();
  _angles.clear();
  _dihedrals.clear();
  _pos = Vec3();
}

// Iterator methods
Bond IXAtom::Begin(AtomBondIter &it) {
  for (it = _bonds.begin(); it != _bonds.end(); ++it) {
    if (!it->expired()) break;
  }
  return (it == _bonds.end()) ? Bond() : it->lock();
}

Bond IXAtom::Next(AtomBondIter &it) {
  for (++it; it != _bonds.end(); ++it) {
    if (!it->expired()) break;
  }
  return (it == _bonds.end()) ? Bond() : it->lock();
}

AtomBondIter IXAtom::BeginBond() { return _bonds.begin(); }
AtomBondIter IXAtom::EndBond() { return _bonds.end(); }

