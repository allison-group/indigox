//
//  pybond.cpp
//  indigox
//
//  Created by Welsh, Ivan on 8/01/18.
//  Copyright © 2018 Allison Group. All rights reserved.
//
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "indigox/python/interface.hpp"

#include "indigox/classes/atom.hpp"
#include "indigox/classes/bond.hpp"
#include "indigox/classes/molecule.hpp"

namespace py = pybind11;

namespace indigox {
  void GeneratePyBond(py::module& m) {
    // Bond class
    py::class_<IXBond, Bond>(m, "Bond")
    // Hidden methods
    .def(py::init<>([](){ return Bond(new IXBond()); }))
    .def("__str__", &IXBond::ToString)
    .def("__repr__", &IXBond::ToString)
    // Getters
    .def("GetIndex", &IXBond::GetIndex)
    .def("GetMolecule", &IXBond::GetMolecule)
    .def("GetOrder", &IXBond::GetOrder)
    .def("GetSourceAtom", &IXBond::GetSourceAtom)
    .def("GetTargetAtom", &IXBond::GetTargetAtom)
    // Setters
    .def("SetIndex", &IXBond::SetIndex)
    .def("SetMolecule", &IXBond::SetMolecule)
    .def("SetOrder", &IXBond::SetOrder)
    .def("SetSourceAtom", &IXBond::SetSourceAtom)
    .def("SetTargetAtom", &IXBond::SetTargetAtom)
    ;
    
    // BondOrder enum
    py::enum_<BondOrder>(m, "BondOrder")
    .value("UNDEFINED_BOND", BondOrder::UNDEFINED_BOND)
    .value("SINGLE_BOND", BondOrder::SINGLE_BOND)
    .value("DOUBLE_BOND", BondOrder::DOUBLE_BOND)
    .value("TRIPLE_BOND", BondOrder::TRIPLE_BOND)
    .value("QUADRUPLE_BOND", BondOrder::QUADRUPLE_BOND)
    .value("AROMATIC_BOND", BondOrder::AROMATIC_BOND)
    .value("ONEANDAHALF_BOND", BondOrder::ONEANDAHALF_BOND)
    .value("TWOANDAHALF_BOND", BondOrder::TWOANDAHALF_BOND)
    ;
  }
}
