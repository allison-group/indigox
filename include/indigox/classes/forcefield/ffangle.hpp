//
//  ffangle.hpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_FORCEFIELD_ANGLE_HPP
#define INDIGOX_FORCEFIELD_ANGLE_HPP

#include "../basedata.hpp"
#include "../../defines.hpp"
#include "../../enumtypes.hpp"

namespace indigox {
    class IXForcefield;
}

namespace indigox {
    namespace ff {
        
        class IXAngleFF : public data::IXBaseData {
            friend class indigox::IXForcefield;
            
        private:
            const IXID m_code;
            const IXANGLEFF_TYPE m_type;
            const float m_equilibirum, m_forceConstant;
            
        private:
            IXAngleFF() = default;
            IXAngleFF(IXID code, float equilibirumValue, float forceConstant, IXANGLEFF_TYPE type);
            
        public:
            inline const IXID GetCode() { return m_code; }
            inline const IXANGLEFF_TYPE GetType() { return m_type; }
            inline const float GetEquilibriumValue() { return m_equilibirum; }
            inline const float GetForceConstant() { return m_forceConstant; }
            
        private:
            void checkValidity();
        };
        
    }  // namespace ff
}  // namespace indigox

#endif /* INDIGOX_FORCEFIELD_ANGLE_HPP */
