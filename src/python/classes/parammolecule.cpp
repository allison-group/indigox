#include <indigox/python/interface.hpp>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>

namespace py = pybind11;

void GeneratePyParameterisedMolecule(pybind11::module& m) {
  using namespace indigox;
  py::return_value_policy Ref = py::return_value_policy::reference;
  
  // ===========================================================================
  // == ParamAtom class bindings ===============================================
  // ===========================================================================
  py::class_<ParamAtom>(m, "ParamAtom")
  .def(py::init<>())
  .def(py::init<const ParamAtom&>())
  .def("MappedWith", &ParamAtom::MappedWith)
  .def("ApplyParameterisation", &ParamAtom::ApplyParameterisation)
  .def("NumSourceAtoms", &ParamAtom::NumSourceAtoms)
  .def("GetAtom", &ParamAtom::GetAtom, Ref)
  .def("MeanCharge", &ParamAtom::MeanCharge)
  .def("MedianCharge", &ParamAtom::MeadianCharge)
  .def("StandardDeviationCharge", &ParamAtom::StandardDeviationCharge)
  .def("GetMostCommonType", &ParamAtom::GetMostCommonType)
  .def("GetMappedTypeCounts", &ParamAtom::GetMappedTypeCounts, Ref)
  .def("GetMappedCharges", &ParamAtom::GetMappedCharges, Ref)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &ParamAtom::operator bool)
  .def("__str__", &outstream_operator<ParamAtom>)
  .def("__repr__", &outstream_operator<ParamAtom>)
  ;
  
  
  // ===========================================================================
  // == ParamBond class bindings ===============================================
  // ===========================================================================
  py::class_<ParamBond>(m, "ParamBond")
  .def(py::init<>())
  .def(py::init<const ParamBond&>())
  .def("MappedWith", &ParamBond::MappedWith)
  .def("ApplyParameterisation", &ParamBond::ApplyParameterisation)
  .def("NumSourceBonds", &ParamBond::NumSourceBonds)
  .def("GetAtoms", &ParamBond::GetAtoms, Ref)
  .def("GetBond", &ParamBond::GetBond, Ref)
  .def("GetMostCommonType", &ParamBond::GetMostCommonType)
  .def("GetMappedTypeCounts", &ParamBond::GetMappedTypeCounts, Ref)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &ParamBond::operator bool)
  .def("__str__", &outstream_operator<ParamBond>)
  .def("__repr__", &outstream_operator<ParamBond>)
  ;
  
  
  // ===========================================================================
  // == ParamAngle class bindings ==============================================
  // ===========================================================================
  py::class_<ParamAngle>(m, "ParamAngle")
  .def(py::init<>())
  .def(py::init<const ParamAngle&>())
  .def("MappedWith", &ParamAngle::MappedWith)
  .def("ApplyParameterisation", &ParamAngle::ApplyParameterisation)
  .def("NumSourceAngles", &ParamAngle::NumSourceAngles)
  .def("GetAtoms", &ParamAngle::GetAtoms, Ref)
  .def("GetAngle", &ParamAngle::GetAngle, Ref)
  .def("GetMostCommonType", &ParamAngle::GetMostCommonType)
  .def("GetMappedTypeCounts", &ParamAngle::GetMappedTypeCounts, Ref)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &ParamAngle::operator bool)
  .def("__str__", &outstream_operator<ParamAngle>)
  .def("__repr__", &outstream_operator<ParamAngle>)
  ;
  
  
  // ===========================================================================
  // == ParamDihedral class bindings ===========================================
  // ===========================================================================
  py::class_<ParamDihedral>(m, "ParamDihedral")
  .def(py::init<>())
  .def(py::init<const ParamDihedral&>())
  .def("MappedWith", &ParamDihedral::MappedWith)
  .def("ApplyParameterisation", &ParamDihedral::ApplyParameterisation)
  .def("NumSourceDihedrals", &ParamDihedral::NumSourceDihedral)
  .def("GetParameterisedAtoms", &ParamDihedral::GetParameterisedAtoms, Ref)
  .def("GetDihedral", &ParamDihedral::GetDihedral, Ref)
  .def("GetMostCommonType", &ParamDihedral::GetMostCommonType)
  .def("GetMappedTypeCounts", &ParamDihedral::GetMappedTypeCounts)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &ParamDihedral::operator bool)
  .def("__str__", &outstream_operator<ParamDihedral>)
  .def("__repr__", &outstream_operator<ParamDihedral>)
  ;
  
  
  // ===========================================================================
  // == ParamMolecule class bindings ===========================================
  // ===========================================================================
  py::class_<ParamMolecule>(m, "ParamMolecule")
  .def(py::init<>())
  .def(py::init<const ParamMolecule&>())
  .def(py::init<Molecule&>())
  .def("ApplyParameterisation", &ParamMolecule::ApplyParameteristion)
  .def("GetAtom", &ParamMolecule::GetAtom)
  .def("GetBond", py::overload_cast<Bond&>(&ParamMolecule::GetBond, py::const_))
  .def("GetBond", py::overload_cast<ParamMolecule::PBond>(&ParamMolecule::GetBond, py::const_))
  .def("GetAngle", py::overload_cast<Angle&>(&ParamMolecule::GetAngle, py::const_))
  .def("GetAngle", py::overload_cast<ParamMolecule::PAngle>(&ParamMolecule::GetAngle, py::const_))
  .def("GetDihedral", py::overload_cast<Dihedral&>(&ParamMolecule::GetDihedral))
  .def("GetDihedral", py::overload_cast<ParamMolecule::PDihedral>(&ParamMolecule::GetDihedral))
  .def("GetAtoms", &ParamMolecule::GetAtoms)
  .def("GetBonds", &ParamMolecule::GetBonds)
  .def("GetAngles", &ParamMolecule::GetAngles)
  .def("GetDihedrals", &ParamMolecule::GetDihedrals)
  ;
  
  
}
