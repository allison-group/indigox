//
//  ffatom.hpp
//  indigox
//
//  Created by Welsh, Ivan on 22/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_FORCEFIELD_ATOM_HPP
#define INDIGOX_FORCEFIELD_ATOM_HPP

#include "../basedata.hpp"
#include "../../defines.hpp"
#include "../../enumtypes.hpp"

namespace indigox {
    class IXForcefield;
}

namespace indigox {
    namespace ff {
        
        class IXAtomFF : public data::IXBaseData {
            friend class indigox::IXForcefield;
            
        private:
            const IXID m_code;
            const IXATOMFF_TYPE m_type;
            const str m_name;
            
            
        private:
            IXAtomFF() = default;
            IXAtomFF(IXATOMFF_TYPE type, IXID code, str name);
            
        public:
            inline const IXID GetCode() { return m_code; }
            inline const IXATOMFF_TYPE GetType() { return m_type; }
            inline const str GetName() { return m_name; }
            
        private:
            void checkValidity();
        };
        
    }  // namespace ff
}  // namespace indigox

#endif /* INDIGOX_FORCEFIELD_ATOM_HPP */
