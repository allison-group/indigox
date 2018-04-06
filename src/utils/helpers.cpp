//
//  helpers.cpp
//  indigox
//
//  Created by Welsh, Ivan on 8/01/18.
//  Copyright © 2018 Allison Group. All rights reserved.
//

#include "indigox/api.hpp"
#include "indigox/classes/atom.hpp"
#include "indigox/classes/bond.hpp"
#include "indigox/classes/molecule.hpp"
#include "indigox/utils/helpers.hpp"

namespace indigox {
  /// @returns a shared pointer to a new Atom instance.
  Atom_p CreateAtom() { return Atom_p(new Atom()); }
  
  /// @returns a shared pointer to a new Bond instance.
  Bond_p CreateBond() { return Bond_p(new Bond()); }
  
  /// @returns a shared pointer to a new Molecule instance.
  Molecule_p CreateMolecule() { return Molecule_p(new Molecule()); }
  
}
