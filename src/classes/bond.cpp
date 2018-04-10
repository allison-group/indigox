/** @file bond.cpp
 *  @brief Bond implementation
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */
#include <sstream>

#include "indigox/classes/atom.hpp"
#include "indigox/classes/bond.hpp"
#include "indigox/utils/counter.hpp"

namespace indigox {
  
  /*
   *  Bond implementation
   */
  
  // Initalisation methods
  IXBond::IXBond() : utils::CountableObject<IXBond>(), mol_() {
    atoms_[0] = _Atom();
    atoms_[1] = _Atom();
    idx_ = GetUniqueID();
  }
  
  /** @param a the source atom for this bond.
   *  @param b the target atom for this bond.
   */
  IXBond::IXBond(Atom a, Atom b) : utils::CountableObject<IXBond>(), mol_() {
    atoms_[0] = std::weak_ptr<IXAtom>(a);
    atoms_[1] = std::weak_ptr<IXAtom>(b);
    idx_ = GetUniqueID();
  }
  
  IXBond::~IXBond() { }
  
  // Data retrival methods
  uid_t IXBond::GetIndex() const { return idx_; }
  Molecule IXBond::GetMolecule() const {
    if (!mol_.expired()) return mol_.lock();
    else return Molecule();
  }
  BondOrder IXBond::GetOrder() const { return order_; }
  Atom IXBond::GetSourceAtom() const {
    if (!atoms_[0].expired()) return atoms_[0].lock();
    else return Atom();
  }
  Atom IXBond::GetTargetAtom() const {
    if (!atoms_[1].expired()) return atoms_[1].lock();
    else return Atom();
  }
  std::string IXBond::ToString() const {
    std::stringstream ss;
    ss << "Bond(";
    if (!atoms_[0].expired()) {
      Atom tmp = atoms_[0].lock();
      ss << tmp->GetIndex() << ":\"" << tmp->GetName() << "\"";
    } else {
      ss << "\"None\"";
    }
    ss << ",";
    if (!atoms_[1].expired()) {
      Atom tmp = atoms_[1].lock();
      ss << tmp->GetIndex() << ":\"" << tmp->GetName() << "\"";
    } else {
      ss << "\"None\"";
    }
    ss << ")";
    return ss.str();
  }
  
  // Data modification methods
  /** @details Index is used for internal bookkeeping and should not be
   *  considered stable.
   *  @param i the index to set.
   */
  void IXBond::SetIndex(uid_t i) { idx_ = i; }
  
  /** @param m the molecule to set. */
  void IXBond::SetMolecule(Molecule m) { mol_ = m; }
  
  /** @param o the bond order to set this bond to. */
  void IXBond::SetOrder(BondOrder o) { order_ = o; }
  
  /** @param a the atom to set this bond's source atom to. */
  void IXBond::SetSourceAtom(Atom a) { atoms_[0] = std::weak_ptr<IXAtom>(a); }
  
  /** @param a the atom to set this bond's target atom to. */
  void IXBond::SetTargetAtom(Atom a) { atoms_[1] = std::weak_ptr<IXAtom>(a); }
  
  void IXBond::Clear() {
    mol_.reset();
    atoms_.fill(std::weak_ptr<IXAtom>());
  }
  
  // Iterator methods
  Atom IXBond::Begin(BondAtomIterator &it) {
    it = atoms_.begin();
    return (it == atoms_.end()) ? Atom() : it->lock();
  }
  
  Atom IXBond::Next(BondAtomIterator &it) {
    ++it;
    return (it == atoms_.end()) ? Atom() : it->lock();
  }
  
  BondAtomIterator IXBond::BeginAtom() { return atoms_.begin(); }
  BondAtomIterator IXBond::EndAtom() { return atoms_.end(); }
  
}

