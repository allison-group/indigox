//
//  ffbond.cpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include "ffbond.hpp"

namespace indigox {
    namespace ff {
        
        IXBondFF::IXBondFF(IXID code, float equilibirumValue, float forceConstant, IXBONDFF_TYPE type)
        : m_code(code), m_equilibirum(equilibirumValue), m_forceConstant(forceConstant), m_type(type) {}
        
        
    }  // namespace ff
}  // namespace indigox