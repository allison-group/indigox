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
    py::class_<IXBond, Bond>(m, "Bond")
    // Only construct from Molecule
    // Hidden methods
    .def("__repr__", &IXBond::ToString)
    // Getters
    .def("GetAromaticity", &IXBond::GetAromaticity)
    .def("GetMolecule", &IXBond::GetMolecule)
    .def("GetOrder", &IXBond::GetOrder)
    .def("GetSourceAtom", &IXBond::GetSourceAtom)
    .def("GetStereochemistry", &IXBond::GetStereochemistry)
    .def("GetTag", &IXBond::GetTag)
    .def("GetTargetAtom", &IXBond::GetTargetAtom)
    .def("NumAngles", &IXBond::NumAngles)
    .def("NumAtoms", &IXBond::NumAtoms)
    .def("NumDihedrals", &IXBond::NumDihedrals)
    // Setters
    .def("SetAromaticity", &IXBond::SetAromaticity)
    .def("SetOrder", &IXBond::SetOrder)
    .def("SetStereochemistry", &IXBond::SetStereochemistry)
    .def("SetTag", &IXBond::SetTag)
    // Other
    .def("SwapSourceTarget", &IXBond::SwapSourceTarget)
    .def("ToString", &IXBond::ToString)
    // Iterators
//    .def("IterateAngles", &IXBond::GetAngleIters)
//    .def("IterateAtoms", &IXBond::GetAtomIters)
//    .def("IterateDihedrals", &IXBond::GetDihedralIters)
    // Pickle support
//    .def(py::pickle(&PickleBond, &UnpickleBond))
    ;
    
    // BondOrder enum
    py::enum_<BondOrder>("BondOrder")
    .value("UNDEFINED", BondOrder::UNDEFINED)
    .value("SINGLE", BondOrder::SINGLE)
    .value("DOUBLE", BondOrder::DOUBLE)
    .value("TRIPLE", BondOrder::TRIPLE)
    .value("QUADRUPLE", BondOrder::QUADRUPLE)
    .value("AROMATIC", BondOrder::AROMATIC)
    .value("ONEANDAHALF", BondOrder::ONEANDAHALF)
    .value("TWOANDAHALF", BondOrder::TWOANDAHALF)
    ;
    
    // BondStereo enum
    py::enum_<BondStereo>("BondStereo")
    .value("UNDEFINED", BondStereo::UNDEFINED)
    .value("NONE", BondStereo::NONE)
    .value("E", BondStereo::E)
    .value("Z", BondStereo::Z)
    ;
  }
}
