#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include <indigox/python/interface.hpp>

#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>

namespace py = pybind11;

void GeneratePyDihedral(pybind11::module& m) {
  using namespace indigox;
  
  auto get_atoms = [](const Dihedral& dhd) {
    stdx::quad<Atom, Atom, Atom, Atom> atms = dhd->GetAtoms();
    return std::make_tuple(atms.first, atms.second, atms.third, atms.fourth);
  };
  
  py::class_<IXDihedral, Dihedral>(m, "Dihedral")
  // Only construct from Molecule
  // Hidden methods
  .def("__repr__", [](Dihedral dhd) {
    std::stringstream ss; ss << dhd; return ss.str();
  })
  // Getters
  .def("GetAtoms", get_atoms)
  .def("GetMolecule", &IXDihedral::GetMolecule)
  .def("GetTag", &IXDihedral::GetTag)
  .def("NumAtoms", &IXDihedral::NumAtoms)
  .def("GetUniqueID", &IXDihedral::GetUniqueID)
  // Setters
  .def("SetTag", &IXDihedral::SetTag)
  // Other
  .def("SwapOrder", &IXDihedral::SwapOrder)
  .def("ToString", &IXDihedral::ToString)
  ;
}

