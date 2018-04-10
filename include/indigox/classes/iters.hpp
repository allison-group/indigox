//
//  iters.hpp
//  indigox
//
//  Created by Ivan Welsh on 7/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "atom.hpp"
#include "bond.hpp"
#include "molecule.hpp"

#ifndef INDIGOX_CLASSES_ITERS_HPP
#define INDIGOX_CLASSES_ITERS_HPP

namespace indigox {
  
  enum IteratorType {
    ATOM_NEIGHBOUR,
    BOND_ATOM,
    MOLECULE_ATOM,
    ATOM_BOND,
    MOLECULE_BOND,
    UNDEFINED
  };
  
  class AtomIterator {
  public:
    AtomIterator();
    AtomIterator(Atom);
    AtomIterator(Bond);
    AtomIterator(Molecule);
    AtomIterator(const AtomIterator&);
    ~AtomIterator();
    
    AtomIterator& operator=(const AtomIterator&);
    operator bool() const;
    // Preincrement
    AtomIterator& operator++();
    // Post increment
    AtomIterator operator++(int);
    // Shared pointer to current atom
    Atom operator->() const;
    Atom operator*() const;
    
  private:
    void ResetToOther(const AtomIterator&);
    
  private:
    IteratorType type_;
    MolAtomIterator itMol_;
    Molecule parentMol_;
    BondAtomIterator itBond_;
    Bond parentBond_;
    AtomBondIter itAtom_;
    Atom parentAtom_;
    
    Atom ptr_;
  };
}

#endif /* INDIGOX_CLASSES_ITERS_HPP */
