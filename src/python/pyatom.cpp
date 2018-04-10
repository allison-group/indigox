//
//  pyatom.cpp
//  indigox
//
//  Created by Welsh, Ivan on 8/01/18.
//  Copyright © 2018 Allison Group. All rights reserved.
//

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "indigox/python/interface.hpp"

#include "indigox/classes/atom.hpp"
#include "indigox/classes/molecule.hpp"
#include "indigox/classes/periodictable.hpp"

namespace py = pybind11;

namespace indigox {
  void GeneratePyAtom(pybind11::module& m) {
    py::class_<IXAtom, Atom>(m, "Atom")
    // Hidden methods
    .def(py::init<>([](){ return Atom(new IXAtom()); }))
    .def("__str__", &IXAtom::ToString)
    .def("__repr__", &IXAtom::ToString)
    // Getters
    .def("GetElement", &IXAtom::GetElement)
    .def("GetFormalCharge", &IXAtom::GetFormalCharge)
    .def("GetIndex", &IXAtom::GetIndex)
    .def("GetMolecule", &IXAtom::GetMolecule)
    .def("GetName", &IXAtom::GetName)
    .def("GetX", &IXAtom::GetX)
    .def("GetY", &IXAtom::GetY)
    .def("GetZ", &IXAtom::GetZ)
    .def("GetUniqueID", &IXAtom::GetUniqueID)
    // Setters
    .def("SetElement", py::overload_cast<Element>(&IXAtom::SetElement))
    .def("SetElement", py::overload_cast<std::string>(&IXAtom::SetElement))
    .def("SetElement", py::overload_cast<unsigned int>(&IXAtom::SetElement))
    .def("SetFormalCharge", &IXAtom::SetFormalCharge)
    .def("SetIndex", &IXAtom::SetIndex)
    .def("SetMolecule", &IXAtom::SetMolecule)
    .def("SetName", &IXAtom::SetName)
    .def("SetPosition", &IXAtom::SetPosition)
    // Pickle support
//    .def(py::pickle(// __getstate__
//                    [](const Atom a) {
//                      return py::make_tuple(a->GetName(),
//                                            a->GetElement()->GetAtomicNumber());
//                    },
//                    // __setstate__
//                    [](py::tuple t) {
//                      if (t.size() != 2)
//                        throw std::runtime_error("Invalid state!");
//                      Atom a = Atom(new IXAtom());
//                      a->SetName(t[0].cast<std::string>());
//                      a->SetElement(t[1].cast<unsigned int>());
//                      return a;
//                    }))
    ;
  }
}
