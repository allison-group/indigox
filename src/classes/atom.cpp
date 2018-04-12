#include <algorithm>
#include <cstdint>
#include <memory>
#include <sstream>

#include "indigox/classes/atom.hpp"
#include "indigox/classes/periodictable.hpp"
#include "indigox/utils/counter.hpp"

using namespace indigox;

IXAtom::IXAtom() : utils::CountableObject<IXAtom>(), _mol(), _elem(), _fc(0),
                   _idx(0), _implicitH(0), _name("ATOM"), _pos({0.0,0.0,0.0}),
                   _partial(0.0), _stereo(ACHIRAL), _aromatic(false) { }

IXAtom::IXAtom(Molecule m) : IXAtom() {
  _mol = m;
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

void IXAtom::Clear() {
  _mol.reset();
  _elem.reset();
  _bonds.clear();
  _angles.clear();
  _dihedrals.clear();
  _pos = Vec3();
}


