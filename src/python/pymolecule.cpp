//
//  pymolecule.cpp
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
#include "indigox/classes/periodictable.hpp"

namespace py = pybind11;
namespace indigox {
  void GeneratePyMolecule(py::module& m) {
    py::class_<IXMolecule, Molecule>(m, "Molecule")
    // Hidden methods
    .def(py::init<>([](){ return Molecule(new IXMolecule()); }))
    // Data retrival
    .def("GetAtomIndex", &IXMolecule::GetAtomIndex)
    .def("GetAtomUniqueID", &IXMolecule::GetAtomUniqueID)
    .def("GetBondIndex", &IXMolecule::GetBondIndex)
    .def("GetBond", &IXMolecule::GetBond)
    .def("GetFormula", &IXMolecule::GetFormula)
    .def("GetName", &IXMolecule::GetName)
    .def("GetTotalCharge", &IXMolecule::GetTotalCharge)
    .def("NumAtoms", &IXMolecule::NumAtoms)
    .def("NumBonds", &IXMolecule::NumBonds)
    // Data modification
    .def("SetName", &IXMolecule::SetName)
    .def("SetTotalCharge", &IXMolecule::SetTotalCharge)
    // Molecule modification
    .def("NewAtom", py::overload_cast<>(&IXMolecule::NewAtom))
    .def("NewAtom", py::overload_cast<Element>(&IXMolecule::NewAtom))
    .def("NewAtom", py::overload_cast<uid_t,Element>(&IXMolecule::NewAtom))
    .def("NewBond", &IXMolecule::NewBond)
    .def("RemoveAtom", &IXMolecule::RemoveAtom)
    .def("RemoveBond", &IXMolecule::RemoveBond)
    // Structural determination
    .def("AssignElectrons", &IXMolecule::AssignElectrons)
    .def("ApplyElectronAssignment", &IXMolecule::ApplyElectronAssignment)
    .def("GetMinimumElectronAssignmentScore", &IXMolecule::GetMinimumElectronAssignmentScore)
    // Iterators
    .def("GetAtoms", [](const Molecule m){
      return py::make_iterator(m->BeginAtom(), m->EndAtom());
    }, py::keep_alive<0, 1>())
    .def("GetBonds", [](const Molecule m){
      return py::make_iterator(m->BeginBond(), m->EndBond());
    }, py::keep_alive<0, 1>())
    ;
  }
}

