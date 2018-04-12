#include <cstdint>
#include <sstream>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "indigox/python/interface.hpp"
#include "indigox/classes/periodictable.hpp"

namespace py = pybind11;
namespace indigox {
  void GeneratePyElement(py::module& m) {
    py::class_<IXElement, Element>(m, "Element")
    // No constructor
    .def("__eq__", py::overload_cast<Element, uint8_t>(&operator==))
    .def("__eq__", py::overload_cast<uint8_t, Element>(&operator==))
    .def("__eq__", py::overload_cast<Element, std::string>(&operator==))
    .def("__eq__", py::overload_cast<std::string, Element>(&operator==))
    .def("__eq__", py::overload_cast<Element, Element>(&operator==))
    .def("__ne__", py::overload_cast<Element, uint8_t>(&operator!=))
    .def("__ne__", py::overload_cast<uint8_t, Element>(&operator!=))
    .def("__ne__", py::overload_cast<Element, std::string>(&operator!=))
    .def("__ne__", py::overload_cast<std::string, Element>(&operator!=))
    .def("__ne__", py::overload_cast<Element, Element>(&operator!=))
    .def("__repr__", [](Element e) {
      std::stringstream ss;
      if (e) ss << "Element(" << e->GetName() << ")";
      return ss.str();
    }) // Data print for developers
    .def("GetAtomicMass", &IXElement::GetAtomicMass)
    .def("GetAtomicNumber", &IXElement::GetAtomicNumber)
    .def("GetAtomicRadius", &IXElement::GetAtomicRadius)
    .def("GetCovalentRadius", &IXElement::GetCovalentRadius)
    .def("GetElectronegativity", &IXElement::GetElectronegativity)
    .def("GetGroup", &IXElement::GetGroup)
    .def("GetHypervalentOctet", &IXElement::GetHypervalentOctet)
    .def("GetName", &IXElement::GetName)
    .def("GetOctet", &IXElement::GetOctet)
    .def("GetPeriod", &IXElement::GetPeriod)
    .def("GetSymbol", &IXElement::GetSymbol)
    .def("GetVanDerWaalsRadius", &IXElement::GetVanDerWaalsRadius)
    .def("GetValenceElectronCount", &IXElement::GetValenceElectronCount)
    .def("ToString", &IXElement::ToString)
    // No pickling support
    ;
  }
}
