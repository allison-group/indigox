#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <indigox/python/interface.hpp>

#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>


namespace py = pybind11;

void GeneratePyAthenaeum(py::module& m) {
  using namespace indigox;
  using namespace indigox::graph;
  py::return_value_policy Ref = py::return_value_policy::reference;
  // ===========================================================================
  // == Enum type bindings =====================================================
  // ===========================================================================
  py::enum_<Fragment::OverlapType>(m, "OverlapType")
  .value("GenericOverlap", Fragment::OverlapType::GenericOverlap)
  ;
  
  // ===========================================================================
  // == Fragment class bindings ================================================
  // ===========================================================================
  py::class_<Fragment>(m, "Fragment")
  .def(py::init<>())
  .def(py::init<MolecularGraph&, std::vector<MGVertex>&, std::vector<MGVertex>&>())
  .def(py::init<const Fragment&>())
  .def("GetGraph", &Fragment::GetGraph, Ref)
  .def("GetFragment", &Fragment::GetFragment)
  .def("Size", &Fragment::Size)
  .def("GetOverlap", &Fragment::GetOverlap)
  .def("IsFragmentVertex", &Fragment::IsFragmentVertex)
  .def("IsOverlapVertex", &Fragment::IsOverlapVertex)
  .def("GetAtoms", &Fragment::GetAtoms)
  .def("GetBonds", &Fragment::GetBonds)
  .def("GetAngles", &Fragment::GetAngles)
  .def("GetDihedrals", &Fragment::GetDihedrals)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  ;
  
  // ===========================================================================
  // == Athenaeum class bindings ===============================================
  // ===========================================================================
  py::class_<Athenaeum, std::unique_ptr<Athenaeum>> PyAthenaeum(m, "Athenaeum");
  PyAthenaeum.def(py::init<Forcefield&>())
  .def(py::init<Forcefield&, uint32_t>())
  .def(py::init<Forcefield&, uint32_t, uint32_t>())
  .def("NumFragments", py::overload_cast<>(&Athenaeum::NumFragments, py::const_))
  .def("NumFragments", py::overload_cast<Molecule&>(&Athenaeum::NumFragments, py::const_))
  .def("GetFragments", py::overload_cast<>(&Athenaeum::GetFragments, py::const_), Ref)
  .def("GetFragments", py::overload_cast<Molecule&>(&Athenaeum::GetFragments, py::const_), Ref)
  .def("HasFragments", &Athenaeum::HasFragments)
  .def("GetForcefield", &Athenaeum::GetForcefield)
  .def("IsSelfConsistent", &Athenaeum::IsSelfConsistent)
  .def("SetSelfConsistent", &Athenaeum::SetSelfConsistent)
  .def("AddFragment", &Athenaeum::AddFragment)
  .def("AddAllFragments", &Athenaeum::AddAllFragments)
  ;
  
  // ===========================================================================
  // == Athenaeum settings bindings ============================================
  // ===========================================================================
  using AtSet = Athenaeum::Settings;
  py::class_<AtSet>(PyAthenaeum, "Settings")
  .def_readwrite_static("AtomLimit", &AtSet::AtomLimit)
  .def_readwrite_static("MinimumFragmentSize", &AtSet::MinimumFragmentSize)
  .def_readwrite_static("MaximumFragmentSize", &AtSet::MaximumFragmentSize)
  .def_readwrite_static("FragmentCycles", &AtSet::FragmentCycles)
  .def_readwrite_static("MaximumCycleSize", &AtSet::MaximumCycleSize)
  .def_readwrite_static("DefaultOverlap", &AtSet::DefaultOverlap)
  .def_readwrite_static("DefaultCycleOverlap", &AtSet::DefaultCycleOverlap)
  ;
}
