//
//  ffdihedral.cpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <stdexcept>

#include "ffdihedral.hpp"

namespace indigox {
    namespace ff {
    
        IXDihedralFF::IXDihedralFF(IXID code, float phase, float forceConstant, int multiplicity, IXDIHEDRALFF_TYPE type)
        : m_code(code), m_equilibirum(phase), m_forceConstant(forceConstant), m_multiplicity(multiplicity), m_type(type) {
            
            checkValidity();
            
        }
        
        void IXDihedralFF::checkValidity() {
            bool allOK = false;
            
            if (m_type == IXDIHEDRALFF_TYPE::PROPER) {
                allOK = (m_equilibirum != fINFINITY) && (m_forceConstant != fINFINITY) && (m_multiplicity != iINFINITY);
                if (!allOK) throw std::invalid_argument("Invalid definition for PROPER dihedral type.");
            } else if (m_type == IXDIHEDRALFF_TYPE::IMPROPER) {
                allOK = (m_equilibirum != fINFINITY) && (m_forceConstant != fINFINITY) && (m_multiplicity == iINFINITY);
                if (!allOK) throw std::invalid_argument("Invalid definition for IMPROPER dihedral type.");
            }
        }
        
    }  // namespace ff
}  // namespace indigox