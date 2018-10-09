#include <cstdint>
#include <set>
#include <sstream>
#include <string>

#include <boost/algorithm/string/join.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <indigox/classes/periodictable.hpp>
#include <indigox/python/interface.hpp>
#include <indigox/python/pyopaquecontainers.hpp>

namespace py = pybind11;

void GeneratePyElement(py::module& m) {
  using namespace indigox;
  py::class_<IXElement, Element>(m, "Element")
  // No constructor
  .def("__eq__", py::overload_cast<const Element&, uint8_t>(&operator==))
  .def("__eq__", py::overload_cast<uint8_t, const Element&>(&operator==))
  .def("__eq__", py::overload_cast<const Element&, std::string>(&operator==))
  .def("__eq__", py::overload_cast<std::string, const Element&>(&operator==))
  .def("__eq__", py::overload_cast<const Element&, const Element&>(&operator==))
  .def("__ne__", py::overload_cast<const Element&, uint8_t>(&operator!=))
  .def("__ne__", py::overload_cast<uint8_t, const Element&>(&operator!=))
  .def("__ne__", py::overload_cast<const Element&, std::string>(&operator!=))
  .def("__ne__", py::overload_cast<std::string, const Element&>(&operator!=))
  .def("__ne__", py::overload_cast<const Element&, const Element&>(&operator!=))
  .def("__repr__", [](const Element& e) {
    std::stringstream ss; ss << e; return ss.str();
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

