//
//  pybond.cpp
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

namespace py = pybind11;

namespace indigox {
  void GeneratePyBond(py::module& m) {
    // Bond class
    py::class_<IXBond, Bond> bond(m, "Bond");
    bond
    // Hidden methods
    .def(py::init<>([](){ return Bond(new IXBond()); }))
    .def(py::init<>([](Atom a, Atom b){ return Bond(new IXBond(a,b)); }))
    .def("__repr__", &IXBond::ToString)
    // Getters
    .def("GetAromaticity", &IXBond::GetAromaticity)
    .def("GetIndex", &IXBond::GetIndex)
    .def("GetMolecule", &IXBond::GetMolecule)
    .def("GetOrder", &IXBond::GetOrder)
    .def("GetSourceAtom", &IXBond::GetSourceAtom)
    .def("GetStereochemistry", &IXBond::GetStereochemistry)
    .def("GetTargetAtom", &IXBond::GetTargetAtom)
    .def("NumAngles", &IXBond::NumAngles)
    .def("NumAtoms", &IXBond::NumAtoms)
    .def("NumDihedrals", &IXBond::NumDihedrals)
    .def("ToString", &IXBond::ToString)
    // Setters
    .def("SetAromaticity", &IXBond::SetAromaticity)
    .def("SetIndex", &IXBond::SetIndex)
    .def("SetMolecule", &IXBond::SetMolecule)
    .def("SetOrder", &IXBond::SetOrder)
    .def("SetSourceAtom", &IXBond::SetSourceAtom)
    .def("SetStereochemistry", &IXBond::SetStereochemistry)
    .def("SetTargetAtom", &IXBond::SetTargetAtom)
    // Modification
//    .def("AddAngle", &IXBond::AddAngle)
//    .def("AddDihedral", &IXBond::AddDihedral)
    .def("Clear", &IXBond::Clear)
//    .def("RemoveAngle", &IXBond::RemoveAngle)
//    .def("RemoveDihedral", &IXBond::RemoveDihedral)
    // Iterators
    .def("IterateAngles", &IXBond::GetAngleIters)
    .def("IterateAtoms", &IXBond::GetAtomIters)
    .def("IterateDihedrals", &IXBond::GetDihedralIters)
    // Pickle support
    .def(py::pickle(&PickleBond, &UnpickleBond))
    ;
    
    // BondOrder enum
    py::enum_<IXBond::Order>(bond, "Order")
    .value("UNDEFINED", IXBond::Order::UNDEFINED)
    .value("SINGLE", IXBond::Order::SINGLE)
    .value("DOUBLE", IXBond::Order::DOUBLE)
    .value("TRIPLE", IXBond::Order::TRIPLE)
    .value("QUADRUPLE", IXBond::Order::QUADRUPLE)
    .value("AROMATIC", IXBond::Order::AROMATIC)
    .value("ONEANDAHALF", IXBond::Order::ONEANDAHALF)
    .value("TWOANDAHALF", IXBond::Order::TWOANDAHALF)
    ;
    
    // BondStereo enum
    py::enum_<IXBond::Stereo>(bond, "Stereo")
    .value("UNDEFINED", IXBond::Stereo::UNDEFINED)
    .value("NONE", IXBond::Stereo::NONE)
    .value("E", IXBond::Stereo::E)
    .value("Z", IXBond::Stereo::Z)
    ;
  }
}
