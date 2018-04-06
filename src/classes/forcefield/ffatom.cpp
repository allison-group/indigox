//
//  ffatom.cpp
//  indigox
//
//  Created by Welsh, Ivan on 22/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include "ffatom.hpp"

namespace indigox {
    namespace ff {
        
        IXAtomFF::IXAtomFF(IXATOMFF_TYPE type, IXID code, str name)
        : m_type(type), m_code(code), m_name(name) {
            
            checkValidity();
        }
        
        void IXAtomFF::checkValidity() {
            
        }
        
    }  // namespace ff
}  // namespace indigox