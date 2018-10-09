#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <indigox/python/interface.hpp>
#include <indigox/python/pyopaquecontainers.hpp>

#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/forcefield.hpp>

namespace py = pybind11;

void GeneratePyAthenaeum(py::module& m) {
  using namespace indigox;
  
  py::class_<IXFragment, Fragment>(m, "Fragment")
  .def("GetGraph", &IXFragment::GetGraph)
  .def("GetFragmentVertices", &IXFragment::GetFragmentVertices)
  ;
  
  py::class_<IXAthenaeum, Athenaeum>(m, "Athenaeum")
  .def(py::init<Forcefield>())
  .def(py::init<Forcefield, uint32_t, uint32_t>())
  .def("AddAllFragments", &IXAthenaeum::AddAllFragments)
  .def("AddFragment", &IXAthenaeum::AddFragment)
  .def("NumFragments", py::overload_cast<>(&IXAthenaeum::NumFragments, py::const_))
  .def("NumFragments", py::overload_cast<const Molecule&>(&IXAthenaeum::NumFragments, py::const_))
  .def("GetFragments", py::overload_cast<>(&IXAthenaeum::GetFragments, py::const_))
  .def("GetFragments", py::overload_cast<const Molecule&>(&IXAthenaeum::GetFragments, py::const_))
  .def("HasFragments", &IXAthenaeum::HasFragments)
  .def("GetForcefield", &IXAthenaeum::GetForcefield)
  .def("IsManualSelfConsistent", &IXAthenaeum::IsManualSelfConsistent)
  .def("SetManual", &IXAthenaeum::SetManual)
  ;
}
