#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>

#include <indigox/classes/periodictable.hpp>

namespace py = pybind11;

void GeneratePyPeriodicTable(py::module& m) {
  using namespace indigox;
  // ===========================================================================
  // == Element class bindings =================================================
  // ===========================================================================
  py::class_<Element>(m, "Element")
  .def(py::init<>())
  .def(py::init<const Element&>())
  .def("GetAtomicMass", &Element::GetAtomicMass)
  .def("GetAtomicNumber", &Element::GetAtomicNumber)
  .def("GetAtomicRadius", &Element::GetAtomicRadius)
  .def("GetCovalentRadius", &Element::GetCovalentRadius)
  .def("GetVanDerWaalsRadius", &Element::GetVanDerWaalsRadius)
  .def("GetName", &Element::GetName)
  .def("GetSymbol", &Element::GetSymbol)
  .def("GetGroup", &Element::GetGroup)
  .def("GetPeriod", &Element::GetPeriod)
  .def("GetValenceElectronCount", &Element::GetValenceElectronCount)
  .def("GetOctet", &Element::GetOctet)
  .def("GetHypervalentOctet", &Element::GetHypervalentOctet)
  .def("GetElectronegativity", &Element::GetElectronegativity)
  .def("__bool__", &Element::operator bool)
  .def(py::self == int32_t())
  .def(py::self == std::string())
  .def(py::self == py::self)
  .def(py::self != int32_t())
  .def(py::self != std::string())
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__str__", &outstream_operator<Element>)
  .def("__repr__", &outstream_operator<Element>)
  ;
  
  // ===========================================================================
  // == PeriodicTable class bindings ===========================================
  // ===========================================================================
  py::class_<PeriodicTable>(m, "PeriodicTable")
  .def("__init__", &GetPeriodicTable, py::return_value_policy::reference)
  .def("GetElement", py::overload_cast<const int32_t>(&PeriodicTable::GetElement, py::const_))
  .def("GetElement", py::overload_cast<const std::string>(&PeriodicTable::GetElement, py::const_))
  .def("__getitem__", py::overload_cast<const int32_t>(&PeriodicTable::operator[], py::const_))
  .def("__getitem__", py::overload_cast<const std::string>(&PeriodicTable::operator[], py::const_))
  .def("GetUndefined", &PeriodicTable::GetUndefined)
  .def("NumElements", &PeriodicTable::NumElements)
  .def("__str__", &outstream_operator<PeriodicTable>)
  .def("__repr__", &outstream_operator<PeriodicTable>)
  ;
  
  // ===========================================================================
  // == Module function bindings ===============================================
  // ===========================================================================
  m.def("GetPeriodicTable", &GetPeriodicTable, py::return_value_policy::reference);
  
}
