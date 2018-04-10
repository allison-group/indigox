//
//  iters.cpp
//  indigox
//
//  Created by Ivan Welsh on 7/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//
#include <iostream>

#include "indigox/classes/iters.hpp"
#include "indigox/classes/molecule.hpp"

using namespace indigox;

AtomIterator::AtomIterator() : type_(UNDEFINED), ptr_() { }

AtomIterator::AtomIterator(Molecule m) : type_(MOLECULE_ATOM) {
  parentMol_ = m;
  ptr_ = parentMol_->Begin(itMol_);
}

AtomIterator::AtomIterator(Bond b) : type_(BOND_ATOM) {
  parentBond_ = b;
  ptr_ = parentBond_->Begin(itBond_);
}

AtomIterator::AtomIterator(Atom a) : type_(ATOM_NEIGHBOUR) {
  parentAtom_ = a;
  Bond tmp = parentAtom_->Begin(itAtom_);
  ptr_ = tmp->GetSourceAtom();
  if (!ptr_ || ptr_->GetUniqueID() != parentAtom_->GetUniqueID()) {
    ptr_ = tmp->GetTargetAtom();
  }
}

AtomIterator::AtomIterator(const AtomIterator& it) {
  ResetToOther(it);
}

AtomIterator::~AtomIterator() { }

AtomIterator& AtomIterator::operator=(const AtomIterator& it) {
  if (this != &it) ResetToOther(it);
  return *this;
}

void AtomIterator::ResetToOther(const AtomIterator& it) {
  type_ = it.type_;
  itMol_ = it.itMol_;
  parentMol_ = it.parentMol_;
  ptr_ = it.ptr_;
}

AtomIterator::operator bool() const { return bool(ptr_); }

AtomIterator& AtomIterator::operator++() {
  switch (type_) {
    case BOND_ATOM:
      ptr_ = parentBond_->Next(itBond_);
      break;
    case MOLECULE_ATOM:
      ptr_ = parentMol_->Next(itMol_);
      break;
    case ATOM_NEIGHBOUR:
      {
        Bond tmp = parentAtom_->Next(itAtom_);
        ptr_ = tmp->GetSourceAtom();
        if (!ptr_ || ptr_->GetUniqueID() != parentAtom_->GetUniqueID()) {
          ptr_ = tmp->GetTargetAtom();
        }
      }
      break;
      
    default:
      break;
  }
  return *this;
}

AtomIterator AtomIterator::operator++(int) {
  AtomIterator tmp(*this);
  operator++();
  return tmp;
}

Atom AtomIterator::operator->() const { return ptr_; }
Atom AtomIterator::operator*() const { return ptr_; }
