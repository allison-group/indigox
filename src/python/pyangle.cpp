#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include <indigox/python/interface.hpp>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/forcefield.hpp>

namespace py = pybind11;

void GeneratePyAngle(pybind11::module& m) {
  using namespace indigox;
  
  auto get_atoms = [](const Angle& ang) {
    stdx::triple<Atom, Atom, Atom> atms = ang->GetAtoms();
    return std::make_tuple(atms.first, atms.second, atms.third);
  };
  
  py::class_<IXAngle, Angle>(m, "Angle")
  // Only construct from Molecule
  // Hidden methods
  .def("__repr__", [](Angle ang) {
    std::stringstream ss; ss << ang; return ss.str();
  })
  // Getters
  .def("GetAtoms", get_atoms)
  .def("GetMolecule", &IXAngle::GetMolecule)
  .def("GetTag", &IXAngle::GetTag)
  .def("NumAtoms", &IXAngle::NumAtoms)
  .def("GetUniqueID", &IXAngle::GetUniqueID)
  .def("GetType", &IXAngle::GetType)
  // Setters
  .def("SetTag", &IXAngle::SetTag)
  .def("SetType", &IXAngle::SetType)
  // Other
  .def("SwapOrder", &IXAngle::SwapOrder)
  .def("ToString", &IXAngle::ToString)
  ;
}
