//
//  forcefield.hpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_FORCEFIELD_HPP
#define INDIGOX_FORCEFIELD_HPP

#include <vector>

#include "basedata.hpp"

#include "../defines.hpp"
#include "../enumtypes.hpp"

namespace indigox {
    namespace ff {
        class IXAngleFF;
        class IXAtomFF;
        class IXBondFF;
        class IXDihedralFF;
    }
    
    class IXForcefield : public data::IXBaseData {
        
    private:
        const IXFORCEFIELDTYPE m_type;
        const str m_version;
        
        std::vector<ff::IXAngleFF*> m_angles;
        std::vector<ff::IXAtomFF*> m_atoms;
        std::vector<ff::IXBondFF*> m_bonds;
        std::vector<ff::IXDihedralFF*> m_dihedrals;
        
        ff::IXANGLEFF_TYPE m_defaultAngleType;
        ff::IXATOMFF_TYPE m_defaultAtomType;
        ff::IXBONDFF_TYPE m_defaultBondType;
        ff::IXDIHEDRALFF_TYPE m_defaultDihedralType;
        
    public:
        IXForcefield();
        IXForcefield(IXFORCEFIELDTYPE type, str version);
        
        ~IXForcefield();
        
        inline void SetDefaultAngleType(ff::IXANGLEFF_TYPE type) { m_defaultAngleType = type; }
        inline void SetDefaultAtomType(ff::IXATOMFF_TYPE type) { m_defaultAtomType = type; }
        inline void SetDefaultBondType(ff::IXBONDFF_TYPE type) { m_defaultBondType = type; }
        inline void SetDefaultDihedralType(ff::IXDIHEDRALFF_TYPE type) { m_defaultDihedralType = type; }
        
        inline ff::IXANGLEFF_TYPE GetDefaultAngleType() const       { return m_defaultAngleType; }
        inline ff::IXATOMFF_TYPE GetDefaultAtomType() const         { return m_defaultAtomType; }
        inline ff::IXBONDFF_TYPE GetDefaultBondType() const         { return m_defaultBondType; }
        inline ff::IXDIHEDRALFF_TYPE GetDefaultDihedralType() const { return m_defaultDihedralType; }
        
        ff::IXAngleFF* NewAngleParameter(IXID code, IXFloat equilibirumValue, IXFloat forceConstant);
        ff::IXAngleFF* NewAngleParameter(ff::IXANGLEFF_TYPE type, IXID code, IXFloat equilibirumValue,
                                         IXFloat forceConstant);
        ff::IXAngleFF* GetAngleParameter(IXID code);
        ff::IXAngleFF* GetAngleParameter(IXID code, ff::IXANGLEFF_TYPE type);
        
        ff::IXAtomFF* NewAtomParameter(IXID code, str name);
        ff::IXAtomFF* NewAtomParameter(ff::IXATOMFF_TYPE type, IXID code, str name);
        ff::IXAtomFF* GetAtomParameter(IXID code);
        ff::IXAtomFF* GetAtomParameter(str name);
        ff::IXAtomFF* GetAtomParameter(IXID code, ff::IXATOMFF_TYPE type);
        ff::IXAtomFF* GetAtomParameter(str name, ff::IXATOMFF_TYPE type);
        
        ff::IXBondFF* NewBondParameter(IXID code, IXFloat equilibirumValue, IXFloat forceConstant);
        ff::IXBondFF* NewBondParameter(ff::IXBONDFF_TYPE type, IXID code, IXFloat equilibirumValue, IXFloat forceConstant);
        ff::IXBondFF* GetBondParameter(IXID code);
        ff::IXBondFF* GetBondParameter(IXID code, ff::IXBONDFF_TYPE type);
        
        // New proper dihedrals
        ff::IXDihedralFF* NewDihedralParameter(IXID code, IXFloat phase, IXFloat forceConstant, int multiplicity);
        ff::IXDihedralFF* NewDihedralParameter(ff::IXDIHEDRALFF_TYPE type, IXID code, IXFloat phase, IXFloat forceConstant, int multiplicity);
        // New improper dihedral
        ff::IXDihedralFF* NewDihedralParameter(IXID code, IXFloat equilibriumValue, IXFloat forceConstant);
        ff::IXDihedralFF* NewDihedralParameter(ff::IXDIHEDRALFF_TYPE type, IXID code, IXFloat equilibriumValue, IXFloat forceConstant);
        ff::IXDihedralFF* GetDihedralParameter(IXID code);
        ff::IXDihedralFF* GetDihedralParameter(IXID code, ff::IXDIHEDRALFF_TYPE type);
        
    private:
        void SetDefaults();
    };
    
}  // namespace indigox

#endif /* INDIGOX_FORCEFIELD_HPP */
