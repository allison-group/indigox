//
//  pyatom.cpp
//  indigox
//
//  Created by Welsh, Ivan on 8/01/18.
//  Copyright © 2018 Allison Group. All rights reserved.
//

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>
#include <indigox/python/pickle.hpp>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>


namespace py = pybind11;

namespace indigox {
  void GeneratePyAtom(pybind11::module& m) {
    py::class_<IXAtom, Atom> atom(m, "Atom");
    atom
    // Constructors
    .def(py::init<>([](){ return Atom(new IXAtom()); }))
    .def(py::init<>([](Molecule m){ return Atom(new IXAtom(m)); }))
    // Hidden methods
    .def("__repr__", &IXAtom::ToString)
    // Getters
    .def("GetAromaticity", &IXAtom::GetAromaticity)
    .def("GetElement", &IXAtom::GetElement)
    .def("GetFormalCharge", &IXAtom::GetFormalCharge)
    .def("GetImplicitCount", &IXAtom::GetImplicitCount)
    .def("GetIndex", &IXAtom::GetIndex)
    .def("GetMolecule", &IXAtom::GetMolecule)
    .def("GetName", &IXAtom::GetName)
    .def("GetPartialCharge", &IXAtom::GetPartialCharge)
    .def("GetStereoChemistry", &IXAtom::GetStereochemistry)
    .def("GetVector", &IXAtom::GetVector)
    .def("GetX", &IXAtom::GetX)
    .def("GetY", &IXAtom::GetY)
    .def("GetZ", &IXAtom::GetZ)
    .def("NumAngles", &IXAtom::NumAngles)
    .def("NumBonds", &IXAtom::NumBonds)
    .def("NumDihedrals", &IXAtom::NumDihedrals)
    .def("GetUniqueID", &IXAtom::GetUniqueID)
    .def("GetCurrentCount", &utils::CountableObject<IXAtom>::GetCurrentCount)
    // Setters
    .def("SetAromaticity", &IXAtom::SetAromaticity)
    .def("SetElement", py::overload_cast<Element>(&IXAtom::SetElement))
    .def("SetElement", py::overload_cast<std::string>(&IXAtom::SetElement))
    .def("SetElement", py::overload_cast<unsigned int>(&IXAtom::SetElement))
    .def("SetFormalCharge", &IXAtom::SetFormalCharge)
    .def("SetImplicitCount", &IXAtom::SetImplicitCount)
    .def("SetIndex", &IXAtom::SetIndex)
    .def("SetMolecule", &IXAtom::SetMolecule)
    .def("SetName", &IXAtom::SetName)
    .def("SetPartialCharge", &IXAtom::SetPartialCharge)
    .def("SetPosition", &IXAtom::SetPosition)
    .def("SetStereochemistry", &IXAtom::SetStereochemistry)
    .def("SetX", &IXAtom::SetX)
    .def("SetY", &IXAtom::SetY)
    .def("SetZ", &IXAtom::SetZ)
    // Modification
//    .def("AddAngle", &IXAtom::AddAngle)
    .def("AddBond", &IXAtom::AddBond)
//    .def("AddDihedral", &IXAtom::AddDihedral)
    .def("Clear", &IXAtom::Clear)
//    .def("RemoveAngle", &IXAtom::RemoveAngle)
    .def("RemoveBond", &IXAtom::RemoveBond)
//    .def("RemoveDihedral", &IXAtom::RemoveDihedral)
    // Other
    .def("ToString", &IXAtom::ToString)
    // Iterators
//    .def("IterateAngles", &IXAtom::GetAngleIters)
    .def("IterateBonds", &IXAtom::GetBondIters)
//    .def("IterateDihedrals", &IXAtom::GetDihedralIters)
    // Pickle support
    .def(py::pickle(&PickleAtom, &UnpickleAtom))
    ;
    
    // AtomStereo enum
    py::enum_<IXAtom::Stereo>(atom, "Stereo")
    .value("UNDEFINED", IXAtom::Stereo::UNDEFINED)
    .value("ACHIRAL", IXAtom::Stereo::ACHIRAL)
    .value("R", IXAtom::Stereo::R)
    .value("S", IXAtom::Stereo::S)
    ;
  }
}
