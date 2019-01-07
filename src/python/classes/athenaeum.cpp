#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/python/interface.hpp>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void GeneratePyAthenaeum(py::module &m) {
  using namespace indigox;
  using namespace indigox::graph;
  py::return_value_policy Ref = py::return_value_policy::reference;
  // ===========================================================================
  // == Enum type bindings =====================================================
  // ===========================================================================
  py::enum_<Fragment::OverlapType>(m, "OverlapType")
      .value("GenericOverlap", Fragment::OverlapType::GenericOverlap);

  // ===========================================================================
  // == Fragment class bindings ================================================
  // ===========================================================================
  py::class_<Fragment>(m, "Fragment")
      .def(py::init<>())
      .def(py::init<const MolecularGraph &, std::vector<MGVertex> &,
                    std::vector<MGVertex> &>())
      .def(py::init<const Molecule &, std::vector<Atom> &,
                    std::vector<Atom> &>())
      .def("GetGraph", &Fragment::GetGraph)
      .def("GetFragment", &Fragment::GetFragment, Ref)
      .def("Size", &Fragment::Size)
      .def("GetOverlap", &Fragment::GetOverlap, Ref)
      .def("IsFragmentVertex", &Fragment::IsFragmentVertex)
      .def("IsOverlapVertex", &Fragment::IsOverlapVertex)
      .def("GetAtoms", &Fragment::GetAtoms, Ref)
      .def("GetBonds", &Fragment::GetBonds, Ref)
      .def("GetAngles", &Fragment::GetAngles, Ref)
      .def("GetDihedrals", &Fragment::GetDihedrals, Ref)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def("__bool__", &Fragment::operator bool);

  // ===========================================================================
  // == Athenaeum class bindings ===============================================
  // ===========================================================================
  py::class_<Athenaeum> PyAthenaeum(m, "Athenaeum");
  PyAthenaeum.def(py::init<>())
      .def(py::init<Forcefield &>())
      .def(py::init<Forcefield &, uint32_t>())
      .def(py::init<Forcefield &, uint32_t, uint32_t>())
      .def("NumFragments",
           py::overload_cast<>(&Athenaeum::NumFragments, py::const_))
      .def("NumFragments", py::overload_cast<const Molecule &>(
                               &Athenaeum::NumFragments, py::const_))
      .def("GetFragments",
           py::overload_cast<>(&Athenaeum::GetFragments, py::const_), Ref)
      .def("GetFragments",
           py::overload_cast<const Molecule &>(&Athenaeum::GetFragments,
                                               py::const_),
           Ref)
      .def("HasFragments", &Athenaeum::HasFragments)
      .def("GetForcefield", &Athenaeum::GetForcefield)
      .def("IsSelfConsistent", &Athenaeum::IsSelfConsistent)
      .def("SetSelfConsistent", &Athenaeum::SetSelfConsistent)
      .def("AddFragment", &Athenaeum::AddFragment)
      .def("AddAllFragments", &Athenaeum::AddAllFragments)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def("__bool__", &Athenaeum::operator bool);

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
      .def_readwrite_static("DefaultCycleOverlap", &AtSet::DefaultCycleOverlap);

  // ===========================================================================
  // == Module level function bindings =========================================
  // ===========================================================================
  m.def("SaveAthenaeum", &SaveAthenaeum);
  m.def("LoadAthenaeum", &LoadAthenaeum);

  // Container bindings
  py::bind_vector<std::vector<Fragment>>(m, "VecFragment");
  py::bind_map<std::map<Molecule, std::vector<Fragment>>>(
      m, "MapMoleculeVecFragment");
}
