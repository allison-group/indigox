#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include <indigox/python/interface.hpp>
#include <indigox/python/pickle.hpp>

#include <indigox/classes/bond.hpp>

namespace py = pybind11;

void GeneratePyBond(py::module& m) {
  using namespace indigox;
  // Bond class
  py::class_<IXBond, Bond>(m, "Bond")
  // Only construct from Molecule
  // Hidden methods
  .def("__repr__", [](Bond bnd) {
    std::stringstream ss; ss << bnd; return ss.str();
  })
  // Getters
  .def("GetAromaticity", &IXBond::GetAromaticity)
  .def("GetAtoms", &IXBond::GetAtoms)
  .def("GetMolecule", &IXBond::GetMolecule)
  .def("GetOrder", &IXBond::GetOrder)
  .def("GetSourceAtom", &IXBond::GetSourceAtom)
  .def("GetStereochemistry", &IXBond::GetStereochemistry)
  .def("GetTag", &IXBond::GetTag)
  .def("GetTargetAtom", &IXBond::GetTargetAtom)
  .def("NumAtoms", &IXBond::NumAtoms)
  .def("GetUniqueID", &IXBond::GetUniqueID)
  // Setters
  .def("SetAromaticity", &IXBond::SetAromaticity)
  .def("SetOrder", &IXBond::SetOrder)
  .def("SetStereochemistry", &IXBond::SetStereochemistry)
  .def("SetTag", &IXBond::SetTag)
  // Other
  .def("SwapSourceTarget", &IXBond::SwapSourceTarget)
  .def("ToString", &IXBond::ToString)
  // Pickle support
  //    .def(py::pickle(&PickleBond, &UnpickleBond))
  ;
  
  // BondOrder enum
  //    py::enum_<BondOrder>("BondOrder")
  //    .value("UNDEFINED", BondOrder::UNDEFINED)
  //    .value("SINGLE", BondOrder::SINGLE)
  //    .value("DOUBLE", BondOrder::DOUBLE)
  //    .value("TRIPLE", BondOrder::TRIPLE)
  //    .value("QUADRUPLE", BondOrder::QUADRUPLE)
  //    .value("AROMATIC", BondOrder::AROMATIC)
  //    .value("ONEANDAHALF", BondOrder::ONEANDAHALF)
  //    .value("TWOANDAHALF", BondOrder::TWOANDAHALF)
  //    ;
  
  // BondStereo enum
  //    py::enum_<BondStereo>("BondStereo")
  //    .value("UNDEFINED", BondStereo::UNDEFINED)
  //    .value("NONE", BondStereo::NONE)
  //    .value("E", BondStereo::E)
  //    .value("Z", BondStereo::Z)
  //    ;
}

