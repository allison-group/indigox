//
//  elimination_ordering.hpp
//  indigox
//
//  Created by Ivan Welsh on 14/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//

#ifndef INDIGOX_ALGORITHM_FORMALBONDS_ELIMINATION_ORDERING_HPP
#define INDIGOX_ALGORITHM_FORMALBONDS_ELIMINATION_ORDERING_HPP

#include <vector>

#include "../../classes/permutablegraph.hpp"

namespace indigox {
  namespace algorithm {
    
    void RandomOrder(PermutableGraph, ElimOrder&);
    void QuickBBOrder(PermutableGraph, ElimOrder&);
    void MinDegreeOrder(PermutableGraph, ElimOrder&);
    void MinAddEdgesOrder(PermutableGraph, ElimOrder&);
    
  }
}

#endif /* INDIGOX_ALGORITHM_FORMALBONDS_ELIMINATION_ORDERING_HPP */
