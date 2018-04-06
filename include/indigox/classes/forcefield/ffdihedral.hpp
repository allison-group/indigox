//
//  ffdihedral.hpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_FORCEFIELD_DIHEDRAL_HPP
#define INDIGOX_FORCEFIELD_DIHEDRAL_HPP

#include "../basedata.hpp"
#include "../../defines.hpp"
#include "../../enumtypes.hpp"

namespace indigox {
    class IXForcefield;
}

namespace indigox {
    namespace ff {
        
        class IXDihedralFF : public data::IXBaseData {
            friend class indigox::IXForcefield;
            
        private:
            const IXID m_code;
            const float m_equilibirum, m_forceConstant;
            const IXDIHEDRALFF_TYPE m_type;
            const int m_multiplicity;
            
        private:
            IXDihedralFF() = default;
            IXDihedralFF(IXID code, float phase, float forceConstant, int multiplicity, IXDIHEDRALFF_TYPE type);
            
        public:
            inline const IXID GetCode() { return m_code; }
            inline const IXDIHEDRALFF_TYPE GetType() { return m_type; }
            inline const float GetEquilibriumValue() { return m_equilibirum; }
            inline const float GetForceConstant() { return m_forceConstant; }
            inline const int GetMultiplicity() { return m_multiplicity; }
            inline const float GetPhaseAngle() { return m_equilibirum; }
            
        private:
            void checkValidity();
        };
        
    }  // namespace ff
}  // namespace indigox

#endif /* INDIGOX_FORCEFIELD_DIHEDRAL_HPP */
