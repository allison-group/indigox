#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "indigox/python/interface.hpp"

#include "indigox/classes/periodictable.hpp"

namespace py = pybind11;
namespace indigox {
  void GeneratePyPeriodicTable(py::module& m) {
    py::class_<IXPeriodicTable, PeriodicTable>(m, "PeriodicTable")
    .def(py::init<>(&IXPeriodicTable::GetInstance), py::return_value_policy::reference)
    .def("GetInstance", &IXPeriodicTable::GetInstance, py::return_value_policy::reference)
    .def("GetElement", py::overload_cast<uint8_t>(&IXPeriodicTable::GetElement))
    .def("GetElement", py::overload_cast<std::string>(&IXPeriodicTable::GetElement))
    .def("__getitem__", py::overload_cast<uint8_t>(&IXPeriodicTable::GetElement))
    .def("__getitem__", py::overload_cast<std::string>(&IXPeriodicTable::GetElement))
    .def("NumElements", &IXPeriodicTable::NumElements)
    ;
  }
}

