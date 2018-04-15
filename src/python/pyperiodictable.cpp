#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>
#include <indigox/python/pickle.hpp>
#include <indigox/classes/periodictable.hpp>

namespace py = pybind11;
namespace indigox {
  void GeneratePyPeriodicTable(py::module& m) {
    py::class_<IXPeriodicTable, PeriodicTable>(m, "PeriodicTable")
    .def(py::init<>(&IXPeriodicTable::GetInstance), py::return_value_policy::reference)
    // allow [] access to elements
    .def("__getitem__", py::overload_cast<const uint8_t>(&IXPeriodicTable::GetElement, py::const_))
    .def("__getitem__", py::overload_cast<const std::string>(&IXPeriodicTable::GetElement, py::const_))
    // allow len(pt) to give the number of elements
    .def("__len__", &IXPeriodicTable::NumElements)
    .def("__repr__", [](PeriodicTable pt) {
      std::stringstream ss;
      if (pt) ss << "PeriodicTable(" << pt->NumElements() << " elements)";
      return ss.str();
    })
    .def("GetInstance", &IXPeriodicTable::GetInstance, py::return_value_policy::reference)
    .def("GetElement", py::overload_cast<const uint8_t>(&IXPeriodicTable::GetElement, py::const_))
    .def("GetElement", py::overload_cast<const std::string>(&IXPeriodicTable::GetElement, py::const_))
    .def("GetUndefinedElement", &IXPeriodicTable::GetUndefinedElement)
    .def("NumElements", &IXPeriodicTable::NumElements)
    .def("ToString", &IXPeriodicTable::ToString)
    // Pickle support
    .def(py::pickle(&PicklePeriodicTable, &UnpicklePeriodicTable))
    ;
  }
}

