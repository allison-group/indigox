#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include <indigox/python/interface.hpp>
#include <indigox/python/pickle.hpp>

#include <indigox/classes/atom.hpp>

namespace py = pybind11;

void GeneratePyAtom(pybind11::module& m) {
  using namespace indigox;
  py::class_<IXAtom, Atom>(m, "Atom")
  // Only construct from Molecule
  // Hidden methods
  .def("__repr__", [](Atom a) {
    std::stringstream ss; ss << a; return ss.str();
  })
  // Getters
  .def("GetAromaticity", &IXAtom::GetAromaticity)
  .def("GetElement", &IXAtom::GetElement)
  .def("GetFormalCharge", &IXAtom::GetFormalCharge)
  .def("GetImplicitCount", &IXAtom::GetImplicitCount)
  .def("GetMolecule", &IXAtom::GetMolecule)
  .def("GetName", &IXAtom::GetName)
  .def("GetPartialCharge", &IXAtom::GetPartialCharge)
  .def("GetStereochemistry", &IXAtom::GetStereochemistry)
  .def("GetTag", &IXAtom::GetTag)
  .def("GetVector", &IXAtom::GetVector)
  .def("GetX", &IXAtom::GetX)
  .def("GetY", &IXAtom::GetY)
  .def("GetZ", &IXAtom::GetZ)
  .def("NumAngles", &IXAtom::NumAngles)
  .def("NumBonds", &IXAtom::NumBonds)
  .def("NumDihedrals", &IXAtom::NumDihedrals)
  .def("GetUniqueID", &IXAtom::GetUniqueID)
  // Setters
  .def("SetAromaticity", &IXAtom::SetAromaticity)
  .def("SetElement", py::overload_cast<Element>(&IXAtom::SetElement))
  .def("SetElement", py::overload_cast<std::string>(&IXAtom::SetElement))
  .def("SetElement", py::overload_cast<unsigned int>(&IXAtom::SetElement))
  .def("SetFormalCharge", &IXAtom::SetFormalCharge)
  .def("SetImplicitCount", &IXAtom::SetImplicitCount)
  .def("SetName", &IXAtom::SetName)
  .def("SetPartialCharge", &IXAtom::SetPartialCharge)
  .def("SetPosition", &IXAtom::SetPosition)
  .def("SetStereochemistry", &IXAtom::SetStereochemistry)
  .def("SetTag", &IXAtom::SetTag)
  .def("SetX", &IXAtom::SetX)
  .def("SetY", &IXAtom::SetY)
  .def("SetZ", &IXAtom::SetZ)
  // Other
  .def("AddImplicitHydrogen", &IXAtom::AddImplicitHydrogen)
  .def("RemoveImplicitHydrogen", &IXAtom::RemoveImplicitHydrogen)
  .def("ToString", &IXAtom::ToString)
  // Pickle support
  //    .def(py::pickle(&PickleAtom, &UnpickleAtom))
  ;
  
  // AtomStereo enum
  //    py::enum_<AtomStereo>("AtomStereo")
  //    .value("UNDEFINED", AtomStereo::UNDEFINED)
  //    .value("ACHIRAL", AtomStereo::ACHIRAL)
  //    .value("R", AtomStereo::R)
  //    .value("S", AtomStereo::S)
  //    ;
}

