//
//  forcefield.cpp
//  indigox
//
//  Created by Welsh, Ivan on 21/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <iostream>
#include <stdexcept>

#include "forcefield.hpp"
#include "forcefield/ffangle.hpp"
#include "forcefield/ffatom.hpp"
#include "forcefield/ffbond.hpp"
#include "forcefield/ffdihedral.hpp"

#include "../defines.hpp"

namespace indigox {

    /************************
     *                      *
     *     Construction     *
     *                      *
     ************************/
    IXForcefield::IXForcefield()
    : m_type(IXFORCEFIELDTYPE::DEFAULT), m_version("54A7") {
        SetDefaults();
    }
    
    // Normal constructor
    IXForcefield::IXForcefield(IXFORCEFIELDTYPE type, str version)
    : m_type(type), m_version(version) {
        SetDefaults();
    }
    
    // Set default new parameter types based on forcefield type.
    void IXForcefield::SetDefaults() {
        using namespace ff;
        if (m_type == IXFORCEFIELDTYPE::DEFAULT) {
            m_defaultAngleType = IXANGLEFF_TYPE::DEFAULT;
            m_defaultAtomType = IXATOMFF_TYPE::DEFAULT;
            m_defaultBondType = IXBONDFF_TYPE::DEFAULT;
            m_defaultDihedralType = IXDIHEDRALFF_TYPE::DEFAULT;
            m_angles.reserve(256);
            m_atoms.reserve(64);
            m_bonds.reserve(256);
            m_dihedrals.reserve(64);
        }
    }
    
    /************************
     *                      *
     *     Destruction      *
     *                      *
     ************************/
    
    IXForcefield::~IXForcefield() {
        for (int i = (int)m_angles.size(); i >= 0; --i) {
            delete m_angles[i];
        }
        for (int i = (int)m_atoms.size(); i >= 0; --i) {
            delete m_atoms[i];
        }
        for (int i = (int)m_bonds.size(); i >= 0; --i) {
            delete m_bonds[i];
        }
        for (int i = (int)m_dihedrals.size(); i >= 0; --i) {
            delete m_dihedrals[i];
        }
    }
    
    
    /************************
     *                      *
     *     Angles stuff     *
     *                      *
     ************************/
    ff::IXAngleFF* IXForcefield::NewAngleParameter(IXID code, IXFloat equilibirumValue, IXFloat forceConstant) {
        return NewAngleParameter(m_defaultAngleType, code, equilibirumValue, forceConstant);
    }
    
    ff::IXAngleFF* IXForcefield::NewAngleParameter(ff::IXANGLEFF_TYPE type, IXID code, IXFloat equilibirumValue,
                                                   IXFloat forceConstant) {
        using namespace ff;
        
        if (GetAngleParameter(code, type) != nullptr) {
            LOGF("Angle ID %d with type %d already exists.\n", code, type);
            return nullptr;
        }
        
        IXAngleFF* newAngle;
        try {
            newAngle = new IXAngleFF(code, equilibirumValue, forceConstant, type);
        } catch (const std::invalid_argument& e) {
            // Replace with proper logging
            LOGF("Invalid argument in angle construction: '%s'\n", e.what());
            return nullptr;
        }
        
        m_angles.push_back(newAngle);
        return newAngle;
    }
    
    ff::IXAngleFF* IXForcefield::GetAngleParameter(IXID code) {
        return GetAngleParameter(code, m_defaultAngleType);
    }
    
    ff::IXAngleFF* IXForcefield::GetAngleParameter(IXID code, ff::IXANGLEFF_TYPE type) {
        for (int i = 0; i < m_angles.size(); ++i) {
            if (m_angles[i]->GetCode() == code && m_angles[i]->GetType() == type) {
                return m_angles[i];
            }
        }
        
        LOGF("Angle ID %d with type %d does not exist.\n", code, type);
        return nullptr;
    }
    
    
    /************************
     *                      *
     *     Atoms stuff      *
     *                      *
     ************************/
    ff::IXAtomFF* IXForcefield::NewAtomParameter(IXID code, str name) {
        return NewAtomParameter(m_defaultAtomType, code, name);
    }
    
    ff::IXAtomFF* IXForcefield::NewAtomParameter(ff::IXATOMFF_TYPE type, IXID code, str name) {
        using namespace ff;
        
        if (GetAtomParameter(code, type) != nullptr) {
            LOGF("Atom ID %d with type %d already exists.\n", code, type);
            return nullptr;
        }
        
        IXAtomFF* newAtom;
        try {
            newAtom = new IXAtomFF(type, code, name);
        } catch (const std::invalid_argument& e) {
            // Replace with proper logging
            LOGF("Invalid argument in atom construction: '%s'\n", e.what());
            return nullptr;
        }
        
        m_atoms.push_back(newAtom);
        return newAtom;
    }
    
    ff::IXAtomFF* IXForcefield::GetAtomParameter(IXID code) {
        return GetAtomParameter(code, m_defaultAtomType);
    }
    
    ff::IXAtomFF* IXForcefield::GetAtomParameter(str name) {
        return GetAtomParameter(name, m_defaultAtomType);
    }
    
    ff::IXAtomFF* IXForcefield::GetAtomParameter(IXID code, ff::IXATOMFF_TYPE type) {
        for (int i = 0; i < m_atoms.size(); ++i) {
            if (m_atoms[i]->GetCode() == code && m_atoms[i]->GetType() == type) {
                return m_atoms[i];
            }
        }
        
        LOGF("Atom ID %d with type %d does not exist.\n", code, type);
        return nullptr;
    }
    
    ff::IXAtomFF* IXForcefield::GetAtomParameter(str name, ff::IXATOMFF_TYPE type) {
        for (int i = 0; i < m_atoms.size(); ++i) {
            if (m_atoms[i]->GetName() == name && m_atoms[i]->GetType() == type) {
                return m_atoms[i];
            }
        }
        
        LOGF("Atom name %s with type %d does not exist.\n", name.c_str(), type);
        return nullptr;
    }
    
    
    /************************
     *                      *
     *     Bonds stuff      *
     *                      *
     ************************/
    ff::IXBondFF* IXForcefield::NewBondParameter(IXID code, IXFloat equilibirumValue, IXFloat forceConstant) {
        return NewBondParameter(m_defaultBondType, code, equilibirumValue, forceConstant);
    }
    
    ff::IXBondFF* IXForcefield::NewBondParameter(ff::IXBONDFF_TYPE type, IXID code, IXFloat equilibirumValue,
                                                 IXFloat forceConstant) {
        using namespace ff;
        
        if (GetBondParameter(code, type) != nullptr) {
            LOGF("Bond ID %d with type %d already exists.\n", code, type);
            return nullptr;
        }
        
        IXBondFF* newBond;
        try {
            newBond = new IXBondFF(code, equilibirumValue, forceConstant, type);
        } catch (const std::invalid_argument& e) {
            // Replace with proper logging
            LOGF("Invalid argument in bond construction: '%s'\n", e.what());
            return nullptr;
        }
        
        m_bonds.push_back(newBond);
        return newBond;
    }
    
    ff::IXBondFF* IXForcefield::GetBondParameter(IXID code) {
        return GetBondParameter(code, m_defaultBondType);
    }
    
    ff::IXBondFF* IXForcefield::GetBondParameter(IXID code, ff::IXBONDFF_TYPE type){
        for (int i = 0; i < m_bonds.size(); ++i) {
            if (m_bonds[i]->GetCode() == code && m_bonds[i]->GetType() == type) {
                return m_bonds[i];
            }
        }
        
        LOGF("Bond ID %d with type %d does not exist.\n", code, type);
        return nullptr;
    }
    
    
    /************************
     *                      *
     *   Dihedrals stuff    *
     *                      *
     ************************/
    ff::IXDihedralFF* IXForcefield::NewDihedralParameter(IXID code, IXFloat phase, IXFloat forceConstant,
                                                         int multiplicity) {
        return NewDihedralParameter(m_defaultDihedralType, code, phase, forceConstant, multiplicity);
    }
    
    ff::IXDihedralFF* IXForcefield::NewDihedralParameter(IXID code, IXFloat phase, IXFloat forceConstant) {
        return NewDihedralParameter(m_defaultDihedralType, code, phase, forceConstant, iINFINITY);
    }
    
    ff::IXDihedralFF* IXForcefield::NewDihedralParameter(ff::IXDIHEDRALFF_TYPE type, IXID code, IXFloat phase,
                                                         IXFloat forceConstant) {
        return NewDihedralParameter(type, code, phase, forceConstant, iINFINITY);
    }
    
    ff::IXDihedralFF* IXForcefield::NewDihedralParameter(ff::IXDIHEDRALFF_TYPE type, IXID code, IXFloat phase,
                                                         IXFloat forceConstant, int multiplicity) {
        using namespace ff;
        
        if (GetDihedralParameter(code, type) != nullptr) {
            LOGF("Dihedral ID %d with type %d already exists.\n", code, type);
            return nullptr;
        }
        
        IXDihedralFF* newDihedral;
        try {
            newDihedral = new IXDihedralFF(code, phase, forceConstant, multiplicity, type);
        } catch (const std::invalid_argument& e) {
            // Replace with proper logging
            LOGF("Invalid argument in dihedral construction: '%s'\n", e.what());
            return nullptr;
        }
        
        m_dihedrals.push_back(newDihedral);
        return newDihedral;
    }
    
    ff::IXDihedralFF* IXForcefield::GetDihedralParameter(IXID code) {
        return GetDihedralParameter(code, m_defaultDihedralType);
    }
    
    ff::IXDihedralFF* IXForcefield::GetDihedralParameter(IXID code, ff::IXDIHEDRALFF_TYPE type){
        for (int i = 0; i < m_dihedrals.size(); ++i) {
            if (m_dihedrals[i]->GetCode() == code && m_dihedrals[i]->GetType() == type) {
                return m_dihedrals[i];
            }
        }
        
        LOGF("Dihedral ID %d with type %d does not exist.\n", code, type);
        return nullptr;
    }
    
}  // namespace indigox
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

