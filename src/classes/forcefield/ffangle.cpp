//
//  ffangle.cpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include "ffangle.hpp"

namespace indigox {
    namespace ff {
        
        IXAngleFF::IXAngleFF(IXID code, float equilibirumValue, float forceConstant, IXANGLEFF_TYPE type)
        : m_code(code), m_equilibirum(equilibirumValue), m_forceConstant(forceConstant), m_type(type) {}
        
    }  // namespace ff
}  // namespace indigox