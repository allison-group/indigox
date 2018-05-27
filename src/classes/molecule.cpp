#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <sstream>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/counter.hpp>
#include <indigox/utils/numerics.hpp>


namespace indigox {
  
  void IXMolecule::SetPropertyModified(Property prop) {
    typedef std::bitset<static_cast<size_t>(Emergent::NUM_EMERGENTS)> es;
    es mask(0);
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
  
  Atom IXMolecule::GetAtomTag(uid_ tag) const {
    auto pos = std::find_if(_atoms.begin(), _atoms.end(),
                       [tag](Atom a) { return a->GetTag() == tag; });
    return (pos == _atoms.end()) ? Atom() : *pos;
  }
  
  Atom IXMolecule::GetAtomID(uid_ id) const {
    auto pos = std::find_if(_atoms.begin(), _atoms.end(),
                            [id](Atom a) { return a->GetUniqueID() == id; });
    return (pos == _atoms.end()) ? Atom() : *pos;
  }
  
  Bond IXMolecule::GetBond(Atom a, Atom b) const {
    IXAtom::AtomBondIter begin, end;
    std::tie(begin, end) = a->GetBondIters();
    auto pos = std::find_if(begin, end, [a,b](_Bond bnd) {
      Bond b_ = bnd.lock();
      if (!b_) return false;
      return ((b_->GetSourceAtom() == a && b_->GetTargetAtom() == b)
              || (b_->GetSourceAtom() == b && b_->GetTargetAtom() == a));
    });
    return (pos == end) ? Bond() : pos->lock();
  }
  
  Bond IXMolecule::GetBondTag(uid_ tag) const {
    auto pos = std::find_if(_bonds.begin(), _bonds.end(),
                            [tag](Bond b) { return  b->GetTag() == tag; });
    return (pos == _bonds.end()) ? Bond() : *pos;
  }
  
  Bond IXMolecule::GetBondID(uid_ id) const {
    auto pos = std::find_if(_bonds.begin(), _bonds.end(),
                            [id](Bond b) { return b->GetUniqueID() == id; });
    return (pos == _bonds.end()) ? Bond() : *pos;
  }
  
  string_ IXMolecule::GetFormula() {
    if (_emerge[static_cast<size_>(Emergent::MOLECULAR_FORMULA)]) {
      std::map<std::string, size_t> e_count;
      for (Atom atm : _atoms) e_count[atm->GetElement()->GetSymbol()]++;
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
  
  size_ IXMolecule::NumAngles() {
    if (_emerge[static_cast<size_>(Emergent::ANGLE_PERCEPTION)]) {
      // DetermineAngles();
      _emerge.reset(static_cast<size_>(Emergent::ANGLE_PERCEPTION));
    }
    return _angles.size();
  }
  
  std::pair<IXMolecule::MolAngleIter, IXMolecule::MolAngleIter>
  IXMolecule::GetAngles() {
    if (_emerge[static_cast<size_>(Emergent::ANGLE_PERCEPTION)]) {
      // DetermineAngles();
      _emerge.reset(static_cast<size_>(Emergent::ANGLE_PERCEPTION));
    }
    return {_angles.cbegin(), _angles.cend()};
  }
  
  size_ IXMolecule::NumDihedrals() {
    if (_emerge[static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION)]) {
      // DetermineDihedrals();
      _emerge.reset(static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION));
    }
    return _dihedrals.size();
  }
  
  std::pair<IXMolecule::MolDihedralIter, IXMolecule::MolDihedralIter>
  IXMolecule::GetDihedrals() {
    if (_emerge[static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION)]) {
      // DetermineDihedrals();
      _emerge.reset(static_cast<size_>(Emergent::DIHEDRAL_PERCEPTION));
    }
    return {_dihedrals.cbegin(), _dihedrals.cend()};
  }
  
  void IXMolecule::SetMolecularCharge(int q) {
    if (q != _q) SetPropertyModified(Property::ELECTRON_COUNT);
    _q = q;
  }
  
  bool IXMolecule::HasAtom(Atom atom) const {
    auto pos = std::find(_atoms.begin(), _atoms.end(), atom);
    return pos != _atoms.end();
  }
  
  bool IXMolecule::HasBond(Bond bond) const {
    auto pos = std::find(_bonds.begin(), _bonds.end(), bond);
    return pos != _bonds.end();
  }
  
  Atom IXMolecule::NewAtom() {
    Atom atom = Atom(new IXAtom(shared_from_this()));
    _atoms.emplace_back(atom);
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
    _bonds.push_back(bond);
    a->AddBond(bond);
    b->AddBond(bond);
    _g->AddEdge(bond);
    SetPropertyModified(Property::CONNECTIVITY);
    return  bond;
  }
  
  bool IXMolecule::RemoveAtom(Atom atom) {
    if (!atom || !HasAtom(atom)) return false;
    
    // Remove all bonds this atom is part of from molecule
    std::vector<Bond> atm_bnds; atm_bnds.reserve(atom->NumBonds());
    if (atom->NumBonds()) SetPropertyModified(Property::CONNECTIVITY);
    auto bnd_itr = atom->GetBondIters();
    for (; bnd_itr.first != bnd_itr.second; ++ bnd_itr.first)
      atm_bnds.emplace_back(bnd_itr.first->lock());
    for (Bond b : atm_bnds) {
      _bonds.erase(std::find(_bonds.begin(), _bonds.end(), b));
      _g->RemoveEdge(_g->GetEdge(b));
    }
    
    // Remove all angles this atom is part of from molecule
    // Remove all dihedrals this atom is part of from molecule
    
    // Remove the atom from the molecule
    _atoms.erase(std::find(_atoms.begin(), _atoms.end(), atom));
    _g->RemoveVertex(_g->GetVertex(atom));
    SetPropertyModified(Property::ATOM_ELEMENTS);
    
    // Clear the atom to invalidate it
    atom->Clear();
    return true;
  }
  
  bool IXMolecule::RemoveBond(Bond bond) {
    if (!bond || !HasBond(bond)) return false;
    
    // Remove the bond from each atom
    bond->GetSourceAtom()->RemoveBond(bond);
    bond->GetTargetAtom()->RemoveBond(bond);
    
    // Remove all angles this bond is aprt of from molecule
    // Remove all dihedrals this bond is part of from molecule
    
    // Remove the bond from the molecule
    _bonds.erase(std::find(_bonds.begin(), _bonds.end(), bond));
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


