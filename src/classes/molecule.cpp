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
  
  Bond IXMolecule::_FindBond(const Atom &a, const Atom &b) const {
    auto pred = [a,b](Bond bnd) {
      std::pair<Atom, Atom> atms = bnd->GetAtoms();
      return ((atms.first == a && atms.second == b)
              || (atms.first == b && atms.second == a));
    };
    MolBondIter pos = std::find_if(_bnds.begin(), _bnds.end(), pred);
    return pos == _bnds.end() ? Bond() : *pos;
  }
  
  Angle IXMolecule::_FindAngle(const Atom &a, const Atom &b, const Atom &c) const {
    stdx::triple<Atom, Atom, Atom> a2c = stdx::make_triple(a, b, c);
    stdx::triple<Atom, Atom, Atom> c2a = stdx::make_triple(c, b, a);
    auto pred = [&a2c, &c2a] (Angle ang) {
      stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
      return (atms == a2c || atms == c2a);
    };
    MolAngleIter pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return pos == _angs.end() ? Angle() : *pos;
  }
  
  Dihedral IXMolecule::_FindDihedral(const Atom &a, const Atom &b,
                                     const Atom &c, const Atom &d) const {
    stdx::quad<Atom, Atom, Atom, Atom> a2d = stdx::make_quad(a, b, c, d);
    stdx::quad<Atom, Atom, Atom, Atom> d2a = stdx::make_quad(d, c, b, a);
    auto pred = [&a2d, &d2a] (Dihedral dhd) {
      stdx::quad<Atom, Atom, Atom, Atom> atms = dhd->GetAtoms();
      return (atms == a2d || atms == d2a);
    };
    MolDihedralIter pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    return pos == _dhds.end() ? Dihedral() : *pos;
  }
  
  Angle IXMolecule::GetAngle(const Atom& a, const Atom&b, const Atom& c) {
    PerceiveAngles();
    return _FindAngle(a, b, c);
  }
  
  Dihedral IXMolecule::GetDihedral(const Atom &a, const Atom &b,
                                   const Atom &c, const Atom &d) {
    PerceiveDihedrals();
    return _FindDihedral(a, b, c, d);
  }
  
  Angle IXMolecule::GetAngleTag(uid_ tag) const {
    auto pred = [tag](Angle ang) { return ang->GetTag() == tag; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return (pos == _angs.end()) ? Angle() : *pos;
  }
  
  Dihedral IXMolecule::GetDihedralTag(uid_ tag) const {
    auto pred = [tag](Dihedral dhd) { return dhd->GetTag() == tag; };
    auto pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    return (pos == _dhds.end()) ? Dihedral() : *pos;
  }
  
  Angle IXMolecule::GetAngleID(uid_ id) const {
    auto pred = [id](Angle ang) { return ang->GetUniqueID() == id; };
    auto pos = std::find_if(_angs.begin(), _angs.end(), pred);
    return (pos == _angs.end()) ? Angle() : *pos;
  }
  
  Dihedral IXMolecule::GetDihedralID(uid_ id) const {
    auto pred = [id](Dihedral dhd) { return dhd->GetUniqueID() == id; };
    auto pos = std::find_if(_dhds.begin(), _dhds.end(), pred);
    return (pos == _dhds.end()) ? Dihedral() : *pos;
  }
  
  size_ IXMolecule::PerceiveAngles() {
    using namespace indigox::graph;
    if (!_emerge[static_cast<size_>(Emergent::ANGLE_PERCEPTION)]) return 0;
    
    // Expected number of angles
    auto sum = [&](size_ current, Atom v) -> size_ {
      size_ degree = v->NumBonds();
      if (degree < 2) return 0;
      return degree * (degree - 1) / 2;
    };
    size_ count = std::accumulate(_atms.begin(), _atms.end(), 0, sum);
    _angs.reserve(count);
    count = 0;
    
    // Adding new angles
    for (MolAtomIter b = _atms.begin(), e = _atms.end(); b != e; ++b) {
      if ((*b)->NumBonds() < 2) continue;
      auto nbrs_it = _g->GetNeighbours(_g->GetVertex(*b));
      std::vector<MGVertex> nbrs(nbrs_it.first, nbrs_it.second);
      for (size_ i = 0; i < nbrs.size() - 1; ++i) {
        Atom a = nbrs[i]->GetAtom();
        for (size_ j = i + 1; j < nbrs.size(); ++j) {
          Atom c = nbrs[j]->GetAtom();
          if (HasAngle(a, *b, c)) continue;
          NewAngle(a, *b, c);
          ++count;
        }
      }
    }
    _emerge.reset(static_cast<size_>(Emergent::ANGLE_PERCEPTION));
    return count;
  }
  
  size_ IXMolecule::PerceiveDihedrals() {
    using namespace indigox::graph;
    if (!_emerge[static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION)]) return 0;
    
    // Expected number of dihedrals
    auto sum = [&](size_ current, Bond b) -> size_ {
      size_ b_degree = b->GetSourceAtom()->NumBonds();
      if (b_degree < 2) return 0;
      size_ c_degree = b->GetTargetAtom()->NumBonds();
      if (c_degree < 2) return 0;
      return (b_degree - 1) * (c_degree - 1);
    };
    size_ count = std::accumulate(_bnds.begin(), _bnds.end(), 0, sum);
    _dhds.reserve(count);
    count = 0;
    
    // Adding new dihedrals
    for (MolBondIter b = _bnds.begin(), e = _bnds.end(); b != e; ++b) {
      Atom B, C;
      std::tie(B,C) = (*b)->GetAtoms();
      if (B->NumBonds() < 2 || C->NumBonds() < 2) continue;
      auto nbrs_it = _g->GetNeighbours(_g->GetVertex(B));
      std::vector<MGVertex> As(nbrs_it.first, nbrs_it.second);
      nbrs_it = _g->GetNeighbours(_g->GetVertex(C));
      std::vector<MGVertex> Ds(nbrs_it.first, nbrs_it.second);
      for (size_ i = 0; i < As.size(); ++i) {
        Atom A = As[i]->GetAtom();
        if (A == C) continue;
        for (size_ j = 0; j < Ds.size(); ++j) {
          Atom D = Ds[j]->GetAtom();
          if (D == B) continue;
          if (HasDihedral(A, B, C, D)) continue;
          NewDihedral(A, B, C, D);
          ++count;
        }
      }
    }
    _emerge.reset(static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION));
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
    return _FindBond(a, b);
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
  
  bool IXMolecule::HasBond(const Atom& a, const Atom& b) const {
    if (!HasAtom(a) || !HasAtom(b)) return false;
    return bool(_FindBond(a, b));
  }
  
  bool IXMolecule::HasAngle(const Atom &a, const Atom &b, const Atom &c) {
    if (!HasAtom(a) || !HasAtom(b) || !HasAtom(c)) return false;
    PerceiveAngles();
    return bool(_FindAngle(a, b, c));
  }
  
  bool IXMolecule::HasDihedral(const Atom &a, const Atom &b,
                               const Atom &c, const Atom &d) {
    if (!HasAtom(a) || !HasAtom(b) || !HasAtom(c) || !HasAtom(d)) return false;
    PerceiveDihedrals();
    return bool(_FindDihedral(a, b, c, d));
  }
  
  Angle IXMolecule::NewAngle(const Atom &a, const Atom &b, const Atom &c) {
    // No need for logic checks as no ability for user to add angles
    _angs.emplace_back(new indigox::IXAngle(a, b, c, shared_from_this()));
    a->AddAngle(_angs.back());
    b->AddAngle(_angs.back());
    c->AddAngle(_angs.back());
    return _angs.back();
  }
  
  Dihedral IXMolecule::NewDihedral(const Atom &a, const Atom &b,
                                   const Atom &c, const Atom &d) {
    // No need for logic checks as no ability for user to add dihedrals
    _dhds.emplace_back(new indigox::IXDihedral(a, b, c, d, shared_from_this()));
    a->AddDihedral(_dhds.back());
    b->AddDihedral(_dhds.back());
    c->AddDihedral(_dhds.back());
    d->AddDihedral(_dhds.back());
    return _dhds.back();
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


