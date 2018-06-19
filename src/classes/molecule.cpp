#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/quad.hpp>
#include <indigox/utils/triple.hpp>


namespace indigox {
  
  void IXMolecule::SetPropertyModified(Property prop) {
    std::bitset<static_cast<size_t>(Emergent::NUM_EMERGENTS)> mask(0);
    switch (prop) {
      case Property::ATOM_ELEMENTS:
        mask.set(static_cast<size_>(Emergent::MOLECULAR_FORMULA));
        mask.set(static_cast<size_>(Emergent::TOPOLOGICAL_BOFC));
        break;
      case Property::CONNECTIVITY:
        mask.set(static_cast<size_>(Emergent::ANGLE_PERCEPTION));
        mask.set(static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION));
        mask.set(static_cast<size_>(Emergent::TOPOLOGICAL_BOFC));
        break;
      case Property::ELECTRON_COUNT:
        mask.set(static_cast<size_>(Emergent::TOPOLOGICAL_BOFC));
        break;
      case Property::NUM_PROPERTIES:
      default:
        break;
    }
    _emerge |= mask;
  }
  
  IXMolecule::IXMolecule() : utils::IXCountableObject<IXMolecule>(), _name(""),
    _q(0), _emerge(0) {
      _emerge.set();
    }
  
  void IXMolecule::Init() {
    _g = graph::MolecularGraph(new graph::IXMolecularGraph(shared_from_this()));
  }
  
  Angle IXMolecule::GetAngle(const Atom& a, const Atom&b, const Atom& c) {
    PerceiveAngles();
    auto pred = [a,b,c](Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return ((atms.first == a && atms.second == b && atms.third == c)
              || (atms.third == a && atms.second == b && atms.first == c));
    };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return (pos == _angs.end()) ? Angle() : *pos;
  }
  
  Angle IXMolecule::GetAngleTag(uid_ tag) const {
    auto pred = [tag](Angle ang) { return ang->GetTag() == tag; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return (pos == _angs.end()) ? Angle() : *pos;
  }
  
  Angle IXMolecule::GetAngleID(uid_ id) const {
    auto pred = [id](Angle ang) { return ang->GetUniqueID() == id; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return (pos == _angs.end()) ? Angle() : *pos;
  }
  
  size_ IXMolecule::PerceiveAngles() {
    if (!_emerge[static_cast<size_>(Emergent::ANGLE_PERCEPTION)]) return 0;
    auto verts = _g->GetVertices();
    auto sum = [&](size_ current, graph::MGVertex v) -> size_ {
      size_ degree = _g->Degree(v);
      if (degree < 2) return 0;
      return degree * (degree - 1) / 2;
    };
    // Expected number of angles
    size_ count = std::accumulate(verts.first, verts.second, 0, sum);
    _angs.reserve(count);
    count = 0;
    
    for (; verts.first != verts.second; ++verts.first) {
      graph::MGVertex b = *verts.first;
      if (_g->Degree(b) < 2) continue;
      auto nbrs_it = _g->GetNeighbours(b);
      std::vector<graph::MGVertex> nbrs(nbrs_it.first, nbrs_it.second);
      for (size_ i = 0; i < nbrs.size() - 1; ++i) {
        graph::MGVertex a = nbrs[i];
        for (size_ j = i + 1; j < nbrs.size(); ++j) {
          graph::MGVertex c = nbrs[j];
          if (HasAngle(a->GetAtom(), b->GetAtom(), c->GetAtom())) continue;
          NewAngle(a->GetAtom(), b->GetAtom(), c->GetAtom());
          ++count;
        }
      }
    }
    _emerge.reset(static_cast<size_>(Emergent::ANGLE_PERCEPTION));
    return count;
  }
  
  Atom IXMolecule::GetAtomTag(uid_ tag) const {
    auto pred = [tag](Atom a) { return a->GetTag() == tag; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? Atom() : *pos;
  }
  
  Atom IXMolecule::GetAtomID(uid_ id) const {
    auto pred = [id](Atom a) { return a->GetUniqueID() == id; };
    auto pos = std::find_if(_atms.begin(), _atms.end(), pred);
    return (pos == _atms.end()) ? Atom() : *pos;
  }
  
  Bond IXMolecule::GetBond(const Atom& a, const Atom& b) const {
    auto pred = [a,b](Bond bnd) {
      std::pair<Atom, Atom> atms = bnd->GetAtoms();
      return ((atms.first == a && atms.second == b)
              || (atms.first == b && atms.second == a));
    };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return (pos == _bnds.end()) ? Bond() : *pos;
  }
  
  Bond IXMolecule::GetBondTag(uid_ tag) const {
    auto pred = [tag](Bond b) { return b->GetTag() == tag; };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return (pos == _bnds.end()) ? Bond() : *pos;
  }
  
  Bond IXMolecule::GetBondID(uid_ id) const {
    auto pred = [id](Bond b) { return b->GetUniqueID() == id; };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return (pos == _bnds.end()) ? Bond() : *pos;
  }
  
  string_ IXMolecule::GetFormula() {
    if (_emerge[static_cast<size_>(Emergent::MOLECULAR_FORMULA)]) {
      std::map<std::string, size_t> e_count;
      for (Atom atm : _atms) e_count[atm->GetElement()->GetSymbol()]++;
      std::stringstream ss;
      if (e_count["C"]) ss << "C";
      if (e_count["C"] > 1) ss << e_count["C"];
      if (e_count["H"]) ss << "H";
      if (e_count["H"] > 1) ss << e_count["H"];
      for (auto& e : e_count) {
        if (e.first != "C" && e.first != "H") {
          ss << e.first;
          if (e.second > 1) ss << e.second;
        }
      }
      _formula_cache = ss.str();
      _emerge.reset(static_cast<size_>(Emergent::MOLECULAR_FORMULA));
    }
    return _formula_cache;
  }
  
//  size_ IXMolecule::NumDihedrals() {
//    if (_emerge[static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION)]) {
//      // DetermineDihedrals();
//      _emerge.reset(static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION));
//    }
//    return _dihedrals.size();
//  }
  
//  std::pair<IXMolecule::MolDihedralIter, IXMolecule::MolDihedralIter>
//  IXMolecule::GetDihedrals() {
//    if (_emerge[static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION)]) {
//      // DetermineDihedrals();
//      _emerge.reset(static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION));
//    }
//    return {_dihedrals.cbegin(), _dihedrals.cend()};
//  }
  
  void IXMolecule::SetMolecularCharge(int q) {
    if (q != _q) SetPropertyModified(Property::ELECTRON_COUNT);
    _q = q;
  }
  
  bool IXMolecule::HasAtom(const Atom& atom) const {
    if (!atom) return false;
    return atom->GetMolecule() == shared_from_this();
  }
  
  bool IXMolecule::HasBond(const Bond& bond) const {
    if (!bond) return false;
    return bond->GetMolecule() == shared_from_this();
  }
  
  bool IXMolecule::HasBond(const Atom& a, const Atom& b) const {
    if (!a || !b) return false;
    auto pred = [a,b](Bond bnd) {
      std::pair<Atom, Atom> atms = bnd->GetAtoms();
      return ((atms.first == a && atms.second == b)
              || (atms.second == a && atms.first == b));
    };
    auto pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return (pos != _bnds.end());
  }
  
  bool IXMolecule::HasAngle(const Angle& angle) const {
    if (!angle) return false;
    return angle->GetMolecule() == shared_from_this();
  }
  
  bool IXMolecule::HasAngle(const Atom &a, const Atom &b, const Atom &c) {
    PerceiveAngles();
    auto pred = [a,b,c] (Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return (b == atms.second && ((a == atms.first && c == atms.third)
                                   || (a == atms.third && c == atms.first)));
    };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return pos != _angs.end();
  }
  
  Angle IXMolecule::NewAngle(const Atom &a, const Atom &b, const Atom &c) {
    // No need for logic checks as no ability for user to add angles
    _angs.emplace_back(new indigox::IXAngle(a, b, c, shared_from_this()));
    a->AddAngle(_angs.back());
    b->AddAngle(_angs.back());
    c->AddAngle(_angs.back());
    return _angs.back();
  }
  
  Atom IXMolecule::NewAtom() {
    Atom atom = Atom(new IXAtom(shared_from_this()));
    _atms.emplace_back(atom);
    _g->AddVertex(atom);
    SetPropertyModified(Property::ATOM_ELEMENTS);
    return atom;
  }
  
  Atom IXMolecule::NewAtom(Element element) {
    Atom atom = NewAtom();
    atom->SetElement(element);
    return atom;
  }
  
  Atom IXMolecule::NewAtom(string_ name) {
    Atom atom = NewAtom();
    atom->SetName(name);
    return atom;
  }
  
  Atom IXMolecule::NewAtom(string_ name, Element element) {
    Atom atom = NewAtom();
    atom->SetName(name);
    atom->SetElement(element);
    return atom;
  }
  
  Bond IXMolecule::NewBond(Atom a, Atom b) {
    if (!a || !b || !HasAtom(a) || !HasAtom(b)) return Bond();
    Bond bond = GetBond(a, b);
    if (bond) return Bond();
    bond = Bond(new indigox::IXBond(a, b, shared_from_this()));
    _bnds.push_back(bond);
    a->AddBond(bond);
    b->AddBond(bond);
    _g->AddEdge(bond);
    SetPropertyModified(Property::CONNECTIVITY);
    return  bond;
  }
  
  bool IXMolecule::RemoveAtom(Atom atom) {
    if (!atom || !HasAtom(atom)) return false;
    
    // Remove all bonds this atom is part of from molecule
    if (atom->NumBonds()) SetPropertyModified(Property::CONNECTIVITY);
    auto bnd_pred = [atom](Bond bnd) {  // Predicate checks if atom in bnd
      std::pair<Atom, Atom> atms = bnd->GetAtoms();
      return !(atms.first == atom || atms.second == atom);
    };
    // Partition so can remove bonds from other atoms
    auto bnd_pos = std::partition(_bnds.begin(), _bnds.end(), bnd_pred);
    auto bnd_pos_erase = bnd_pos;
    for (auto it = _bnds.end(); bnd_pos != it; ++bnd_pos) {
      std::pair<Atom, Atom> atms = (*bnd_pos)->GetAtoms();
      if (atms.first == atom) atms.second->RemoveBond(*bnd_pos);
      else atms.first->RemoveBond(*bnd_pos);
    }
    _bnds.erase(bnd_pos_erase, _bnds.end());
    
    // Remove all angles this atom is part of from molecule
    auto ang_pred = [atom](Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return !(atms.first == atom || atms.second == atom || atms.third == atom);
    };
    auto ang_pos = std::partition(_angs.begin(), _angs.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = _angs.end(); ang_pos != it; ++ang_pos) {
      stdx::triple<Atom, Atom, Atom> atms = (*ang_pos)->GetAtoms();
      if (atms.first == atom) {
        atms.second->RemoveAngle(*ang_pos);
        atms.third->RemoveAngle(*ang_pos);
      } else if (atms.second == atom) {
        atms.first->RemoveAngle(*ang_pos);
        atms.third->RemoveAngle(*ang_pos);
      } else {
        atms.second->RemoveAngle(*ang_pos);
        atms.third->RemoveAngle(*ang_pos);
      }
    }
    _angs.erase(ang_pos_erase, _angs.end());
    
    // Remove all dihedrals this atom is part of from molecule
    auto dhd_pred = [atom](Dihedral dhd) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = dhd->GetAtoms();
      return !(atms.first == atom || atms.second == atom
               || atms.third == atom || atms.fourth == atom);
    };
    auto dhd_pos = std::partition(_dhds.begin(), _dhds.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = _dhds.end(); dhd_pos != it; ++dhd_pos) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = (*dhd_pos)->GetAtoms();
      if (atms.first == atom) {
        atms.second->RemoveDihedral(*dhd_pos);
        atms.third->RemoveDihedral(*dhd_pos);
        atms.fourth->RemoveDihedral(*dhd_pos);
      } else if (atms.second == atom) {
        atms.first->RemoveDihedral(*dhd_pos);
        atms.third->RemoveDihedral(*dhd_pos);
        atms.fourth->RemoveDihedral(*dhd_pos);
      } else if (atms.third == atom) {
        atms.second->RemoveDihedral(*dhd_pos);
        atms.first->RemoveDihedral(*dhd_pos);
        atms.fourth->RemoveDihedral(*dhd_pos);
      } else {
        atms.second->RemoveDihedral(*dhd_pos);
        atms.third->RemoveDihedral(*dhd_pos);
        atms.first->RemoveDihedral(*dhd_pos);
      }
    }
    _dhds.erase(dhd_pos_erase, _dhds.end());
    
    // Remove the atom from the molecule
    _atms.erase(std::find(_atms.begin(), _atms.end(), atom));
    _g->RemoveVertex(_g->GetVertex(atom));
    SetPropertyModified(Property::ATOM_ELEMENTS);
    
    // Clear the atom to invalidate it
    atom->Clear();
    return true;
  }
  
  bool IXMolecule::RemoveBond(Bond bond) {
    if (!bond || !HasBond(bond)) return false;
    Atom a, b;
    std::tie(a, b) = bond->GetAtoms();
    // Remove the bond from each atom
    a->RemoveBond(bond);
    b->RemoveBond(bond);
    
    // Remove all angles this bond is part of from molecule
    auto ang_pred = [a,b] (Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return ((a == atms.second && (b == atms.first || b == atms.third))
              || (b == atms.second && (a == atms.first || a == atms.third)));
    };
    auto ang_pos = std::partition(_angs.begin(), _angs.end(), ang_pred);
    auto ang_pos_erase = ang_pos;
    for (auto it = _angs.end(); it != ang_pos; ++ang_pos) {
      stdx::triple<Atom, Atom, Atom> atms = (*ang_pos)->GetAtoms();
      atms.first->RemoveAngle(*ang_pos);
      atms.second->RemoveAngle(*ang_pos);
      atms.third->RemoveAngle(*ang_pos);
    }
    _angs.erase(ang_pos_erase, _angs.end());
    
    // Remove all dihedrals this bond is part of from molecule
    auto dhd_pred = [a,b] (Dihedral dhd) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = dhd->GetAtoms();
      return ((a == atms.second && (b == atms.first || b == atms.third))
              || (a == atms.third && (b == atms.second || b == atms.fourth))
              || (b == atms.second && (a == atms.first || a == atms.third))
              || (b == atms.third && (a == atms.second || a == atms.fourth)));
    };
    auto dhd_pos = std::partition(_dhds.begin(), _dhds.end(), dhd_pred);
    auto dhd_pos_erase = dhd_pos;
    for (auto it = _dhds.end(); dhd_pos != it; ++dhd_pos) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = (*dhd_pos)->GetAtoms();
      atms.first->RemoveDihedral(*dhd_pos);
      atms.second->RemoveDihedral(*dhd_pos);
      atms.third->RemoveDihedral(*dhd_pos);
      atms.fourth->RemoveDihedral(*dhd_pos);
    }
    _dhds.erase(dhd_pos_erase, _dhds.end());
    
    // Remove the bond from the molecule
    _bnds.erase(std::find(_bnds.begin(), _bnds.end(), bond));
    _g->RemoveEdge(_g->GetEdge(bond));
    SetPropertyModified(Property::CONNECTIVITY);
    
    // Clear the bond to invalidate it
    bond->Clear();
    return true;
  }
}


//using namespace indigox;




/*
 * Molecule implementation
 */

//size_t IXMolecule::AssignElectrons() {
//  // check is single component (use graph)
//  /// @todo
//  if (graph_->NumConnectedComponents() != 1) {
//    std::cerr << "Electrons can only be assigned when molecule has only one component." << std::endl;
//    return 0;
//  }
//
//  // check all elements valid.
//  typedef Options::AssignElectrons opt_;
//  for (Atom atm : _atoms) {
//    if (!atm->GetElement()) {
//      std::cerr << "Not all atoms have elements assigned." << std::endl;
//      return 0;
//    }
//    if (opt_::ALLOWED_ELEMENTS.find(atm->GetElement()->GetSymbol())
//        == opt_::ALLOWED_ELEMENTS.end()) {
//      std::cerr << atm->GetElement() << " is not an allowed element for the algorithm." << std::endl;
//      return 0;
//    }
//  }
//
//  // Set the ElectronOpt graph
//  elnopt_->SetMolecularGraph(graph_);
//  modified_ = false;
//  size_t num = (size_t)elnopt_->Run();
//  if (num) ApplyElectronAssignment(0);
//  for (auto& tmp : idx_to_atom_) {
//    Atom_p atom = tmp.second.lock();
//    if (atom->GetFormalCharge() != 0) {
//      std::cout << tmp.first << " : " << atom->GetFormalCharge() << std::endl;
//    }
//  }
//
//  for (auto& tmp : idx_to_bond_) {
//    Bond_p bond = tmp.second.lock();
//    if (bond->GetOrder() != Order::SINGLE_BOND) {
//      uid_t a = bond->GetSourceAtom()->GetIndex();
//      uid_t b = bond->GetTargetAtom()->GetIndex();
//      if (a <= b) std::cout << a << "-" << b << " : " << bond->GetOrder() << std::endl;
//      else std::cout << b << "-" << a << " : " << bond->GetOrder() << std::endl;
//    }
//  }
//  return num;
//}
//
//FCSCORE IXMolecule::GetMinimumElectronAssignmentScore() {
//  if (modified_) {
//    std::cerr << "Molecule has been modified since assignment. Please re-assign." << std::endl;
//    return Options::AssignElectrons::INF;
//  }
//  return elnopt_->GetMinimisedEnergy();
//}
//
//bool IXMolecule::ApplyElectronAssignment(size_t idx) {
//  if (modified_) {
//    std::cerr << "Molecule has been modified since assignment. Please re-assign." << std::endl;
//    return false;
//  }
//  return elnopt_->ApplyElectronAssigment(idx);
//}
//
//// Iterator methods
//Atom IXMolecule::Begin(MolAtomIterator &it) {
//  it = _atoms.begin();
//  return (it == _atoms.end()) ? Atom() : *it;
//}
//
//Atom IXMolecule::Next(MolAtomIterator &it) {
//  ++it;
//  return (it == _atoms.end()) ? Atom() : *it;
//}
//
//MolAtomIterator IXMolecule::BeginAtom() { return _atoms.begin(); }
//MolAtomIterator IXMolecule::EndAtom() { return _atoms.end(); }
//MolBondIterator IXMolecule::BeginBond() { return _bonds.begin(); }
//MolBondIterator IXMolecule::EndBond() { return _bonds.end(); }


