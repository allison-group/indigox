#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>
#include <indigox/python/interface.hpp>

namespace py = pybind11;

// PYBIND11_MAKE_OPAQUE(indigox::ParamAtom::TypeCounts);
// PYBIND11_MAKE_OPAQUE(indigox::ParamBond::TypeCounts);
// PYBIND11_MAKE_OPAQUE(indigox::ParamAngle::TypeCounts);
// PYBIND11_MAKE_OPAQUE(indigox::ParamDihedral::TypeGroup);
// PYBIND11_MAKE_OPAQUE(indigox::ParamDihedral::TypeCounts);

void GeneratePyParameterisedMolecule(pybind11::module &m) {
  using namespace indigox;
  py::return_value_policy Ref = py::return_value_policy::reference;

  // ===========================================================================
  // == ParamAtom class bindings ===============================================
  // ===========================================================================
  py::class_<ParamAtom>(m, "ParamAtom")
      .def(py::init<>())
      .def("MappedWith", &ParamAtom::MappedWith)
      .def("ApplyParameterisation", &ParamAtom::ApplyParameterisation)
      .def("NumSourceAtoms", &ParamAtom::NumSourceAtoms)
      .def("GetAtom", &ParamAtom::GetAtom)
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
      .def("__repr__", &outstream_operator<ParamAtom>);

  // ===========================================================================
  // == ParamBond class bindings ===============================================
  // ===========================================================================
  py::class_<ParamBond>(m, "ParamBond")
      .def(py::init<>())
      .def("MappedWith", &ParamBond::MappedWith)
      .def("ApplyParameterisation", &ParamBond::ApplyParameterisation)
      .def("NumSourceBonds", &ParamBond::NumSourceBonds)
      .def("GetAtoms", &ParamBond::GetAtoms)
      .def("GetBond", &ParamBond::GetBond)
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
      .def("__repr__", &outstream_operator<ParamBond>);

  // ===========================================================================
  // == ParamAngle class bindings ==============================================
  // ===========================================================================
  py::class_<ParamAngle>(m, "ParamAngle")
      .def(py::init<>())
      .def(py::init<const ParamAngle &>())
      .def("MappedWith", &ParamAngle::MappedWith)
      .def("ApplyParameterisation", &ParamAngle::ApplyParameterisation)
      .def("NumSourceAngles", &ParamAngle::NumSourceAngles)
      .def("GetAtoms", &ParamAngle::GetAtoms)
      .def("GetAngle", &ParamAngle::GetAngle)
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
      .def("__repr__", &outstream_operator<ParamAngle>);

  // ===========================================================================
  // == ParamDihedral class bindings ===========================================
  // ===========================================================================
  py::class_<ParamDihedral>(m, "ParamDihedral")
      .def(py::init<>())
      .def(py::init<const ParamDihedral &>())
      .def("MappedWith", &ParamDihedral::MappedWith)
      .def("ApplyParameterisation", &ParamDihedral::ApplyParameterisation)
      .def("NumSourceDihedrals", &ParamDihedral::NumSourceDihedral)
      .def("GetAtoms", &ParamDihedral::GetAtoms)
      .def("GetDihedral", &ParamDihedral::GetDihedral)
      .def("GetMostCommonType", &ParamDihedral::GetMostCommonType)
      .def("GetMappedTypeCounts", &ParamDihedral::GetMappedTypeCounts, Ref)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def("__bool__", &ParamDihedral::operator bool)
      .def("__str__", &outstream_operator<ParamDihedral>)
      .def("__repr__", &outstream_operator<ParamDihedral>);

  // ===========================================================================
  // == ParamMolecule class bindings ===========================================
  // ===========================================================================
  py::class_<ParamMolecule>(m, "ParamMolecule")
      .def(py::init<>())
      .def(py::init<Molecule &>())
      .def("ApplyParameterisation", &ParamMolecule::ApplyParameteristion)
      .def("GetAtom", &ParamMolecule::GetAtom)
      .def("GetBond",
           py::overload_cast<const Bond &>(&ParamMolecule::GetBond, py::const_))
      .def("GetBond", py::overload_cast<const Atom &, const Atom &>(
                          &ParamMolecule::GetBond, py::const_))
      .def("GetAngle", py::overload_cast<const Angle &>(
                           &ParamMolecule::GetAngle, py::const_))
      .def("GetAngle",
           py::overload_cast<const Atom &, const Atom &, const Atom &>(
               &ParamMolecule::GetAngle, py::const_))
      .def("GetDihedral",
           py::overload_cast<const Dihedral &>(&ParamMolecule::GetDihedral))
      .def("GetDihedral",
           py::overload_cast<const Atom &, const Atom &, const Atom &,
                             const Atom &>(&ParamMolecule::GetDihedral))
      .def("GetAtoms", &ParamMolecule::GetAtoms, Ref)
      .def("GetBonds", &ParamMolecule::GetBonds, Ref)
      .def("GetAngles", &ParamMolecule::GetAngles, Ref)
      .def("GetDihedrals", &ParamMolecule::GetDihedrals, Ref);

  // container bindings
  py::bind_vector<std::vector<ParamAtom>>(m, "VecParamAtom");
  py::bind_vector<std::vector<ParamBond>>(m, "VecParamBond");
  py::bind_vector<std::vector<ParamAngle>>(m, "VecParamAngle");
  py::bind_vector<std::vector<ParamDihedral>>(m, "VecParamDihedral");

  py::bind_map<eastl::vector_map<indigox::FFAtom, size_t>>(m, "MapFFAtomUInt");
  py::bind_map<eastl::vector_map<indigox::FFBond, size_t>>(m, "MapFFBondUInt");
  py::bind_map<eastl::vector_map<indigox::FFAngle, size_t>>(m,
                                                            "MapFFAngleUInt");
  py::bind_vector<std::vector<FFDihedral>>(m, "VecFFDihedral");
  py::bind_map<eastl::vector_map<std::vector<indigox::FFDihedral>, size_t>>(
      m, "MapVecFFDihedralUInt");
}
