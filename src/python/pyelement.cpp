#include <cstdint>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "indigox/python/interface.hpp"
#include "indigox/classes/periodictable.hpp"

namespace py = pybind11;
namespace indigox {
  void GeneratePyElement(py::module& m) {
    py::class_<IXElement, Element>(m, "Element")
    .def("__eq__", py::overload_cast<Element, uint8_t>(&operator==))
    .def("__eq__", py::overload_cast<Element, std::string>(&operator==))
    .def("__eq__", py::overload_cast<Element, Element>(&operator==))
    .def("__ne__", py::overload_cast<Element, uint8_t>(&operator!=))
    .def("__ne__", py::overload_cast<Element, std::string>(&operator!=))
    .def("__ne__", py::overload_cast<Element, Element>(&operator!=))
    .def("__str__", &IXElement::ToString)
    .def("__repr__", &IXElement::ToString)
    .def("GetAtomicMass", &IXElement::GetAtomicMass)
    .def("GetAtomicNumber", &IXElement::GetAtomicNumber)
    .def("GetAtomicRadius", &IXElement::GetAtomicRadius)
    .def("GetCovalentRadius", &IXElement::GetCovalentRadius)
    .def("GetVanDerWaalsRadius", &IXElement::GetVanDerWaalsRadius)
    .def("GetName", &IXElement::GetName)
    .def("GetSymbol", &IXElement::GetSymbol)
    .def("GetGroup", &IXElement::GetGroup)
    .def("GetPeriod", &IXElement::GetPeriod)
    .def("GetValenceElectronCount", &IXElement::GetValenceElectronCount)
    .def("GetOctet", &IXElement::GetOctet)
    .def("GetHypervalentOctet", &IXElement::GetHypervalentOctet)
    .def("GetElectronegativity", &IXElement::GetElectronegativity)
    ;
  }
}
