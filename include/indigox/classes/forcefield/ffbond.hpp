//
//  ffbond.hpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_FORCEFIELD_BOND_HPP
#define INDIGOX_FORCEFIELD_BOND_HPP

#include "../basedata.hpp"
#include "../../defines.hpp"
#include "../../enumtypes.hpp"

namespace indigox {
    class IXForcefield;
}

namespace indigox {
    namespace ff {
        
        class IXBondFF : public data::IXBaseData {
            friend class indigox::IXForcefield;
            
        private:
            const IXID m_code;
            const IXBONDFF_TYPE m_type;
            const float m_equilibirum, m_forceConstant;
            
            
        private:
            IXBondFF() = default;
            IXBondFF(IXID code, float equilibirumValue, float forceConstant, IXBONDFF_TYPE type);
            
        public:
            inline const IXID GetCode() { return m_code; }
            inline const IXBONDFF_TYPE GetType() { return m_type; }
            inline const float GetEquilibriumValue() { return m_equilibirum; }
            inline const float GetForceConstant() { return m_forceConstant; }
            
        private:
            void checkValidity();
            
        };
        
    }  // namespace ff
}  // namespace indigox

#endif /* INDIGOX_FORCEFIELD_BOND_HPP */
