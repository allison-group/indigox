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
  py::class_<Element>(m, "Element")
  // No constructor
  .def(py::self == py::self)
  .def(py::self == uint8_t())
  .def(py::self == std::string())
  .def(py::self != py::self)
  .def(py::self != uint8_t())
  .def(py::self != std::string())
  .def("__str__", &Element::ToString)
  .def("__repr__", &Element::ToString)
  .def("GetAtomicMass", &Element::GetAtomicMass)
  .def("GetAtomicNumber", &Element::GetAtomicNumber)
  .def("GetAtomicRadius", &Element::GetAtomicRadius)
  .def("GetCovalentRadius", &Element::GetCovalentRadius)
  .def("GetElectronegativity", &Element::GetElectronegativity)
  .def("GetGroup", &Element::GetGroup)
  .def("GetHypervalentOctet", &Element::GetHypervalentOctet)
  .def("GetName", &Element::GetName)
  .def("GetOctet", &Element::GetOctet)
  .def("GetPeriod", &Element::GetPeriod)
  .def("GetSymbol", &Element::GetSymbol)
  .def("GetVanDerWaalsRadius", &Element::GetVanDerWaalsRadius)
  .def("GetValenceElectronCount", &Element::GetValenceElectronCount)
  .def("ToString", &Element::ToString)
  // No pickling support
  ;
}

