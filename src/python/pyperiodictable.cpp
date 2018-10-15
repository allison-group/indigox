#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>
#include <indigox/python/pickle.hpp>

#include <indigox/classes/periodictable.hpp>

namespace py = pybind11;

void GeneratePyPeriodicTable(py::module& m) {
  using namespace indigox;
  
  m.def("GetPeriodicTable", &GetPeriodicTable);
  
  py::class_<PeriodicTable>(m, "PeriodicTable")
  // No constructor
  // allow [] access to elements
  .def("__getitem__", py::overload_cast<const uint8_t>(&PeriodicTable::GetElement, py::const_))
  .def("__getitem__", py::overload_cast<const std::string>(&PeriodicTable::GetElement, py::const_))
  // allow len(pt) to give the number of elements
  .def("__len__", &PeriodicTable::NumElements)
  .def("__repr__", [](const PeriodicTable& pt) { return pt.ToString(); })
  .def("GetElement", py::overload_cast<const uint8_t>(&PeriodicTable::GetElement, py::const_))
  .def("GetElement", py::overload_cast<const std::string>(&PeriodicTable::GetElement, py::const_))
  .def("GetUndefined", &PeriodicTable::GetUndefined)
  .def("NumElements", &PeriodicTable::NumElements)
  .def("ToString", &PeriodicTable::ToString)
  // Pickle support
//  .def(py::pickle(&PicklePeriodicTable, &UnpicklePeriodicTable))
  ;
}
