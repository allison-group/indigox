#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>

#include <indigox/classes/forcefield.hpp>

namespace py = pybind11;

void GeneratePyForcefield(py::module& m) {
  using namespace indigox;
  // ===========================================================================
  // == Enum type bindings =====================================================
  // ===========================================================================
  py::enum_<FFFamily>(m, "FFFamily")
  .value("Empty", FFFamily::Empty)
  .value("GROMOS", FFFamily::GROMOS)
  .value("CHARMM", FFFamily::CHARMM)
  .value("AMBER", FFFamily::AMBER)
  .value("Other", FFFamily::Other)
  ;
  
  py::enum_<BondType>(m, "BondType")
  .value("Empty", BondType::Empty)
  .value("Harmonic", BondType::Harmonic)
  .value("Quartic", BondType::Quartic)
  .value("Morse", BondType::Morse)
  .value("Cubic", BondType::Cubic)
  .value("FENE", BondType::FENE)
  ;
  
  py::enum_<AngleType>(m, "AngleType")
  .value("Empty", AngleType::Empty)
  .value("Harmonic", AngleType::Harmonic)
  .value("CosineHarmonic", AngleType::CosineHarmonic)
  .value("UreyBradley", AngleType::UreyBradley)
  .value("Quartic", AngleType::Quartic)
  ;
  
  py::enum_<DihedralType>(m, "DihedralType")
  .value("Empty", DihedralType::Empty)
  .value("Proper", DihedralType::Proper)
  .value("Improper", DihedralType::Improper)
  .value("RyckaertBellemans", DihedralType::RyckaertBellemans)
  .value("PeriodicImproper", DihedralType::PeriodicImproper)
  .value("Fourier", DihedralType::Fourier)
  .value("Restricted", DihedralType::Restricted)
  ;
  
  // ===========================================================================
  // == FFAtom class bindings ==================================================
  // ===========================================================================
  py::class_<FFAtom>(m, "FFAtom")
  .def(py::init<>())
  .def("GetID", &FFAtom::GetID)
  .def("GetName", &FFAtom::GetName)
  .def("GetForcefield", &FFAtom::GetForcefield)
  .def("GetElement", &FFAtom::GetElement)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &FFAtom::operator bool)
  .def("__str__", &outstream_operator<FFAtom>)
  .def("__repr__", &outstream_operator<FFAtom>)
  ;
  
  // ===========================================================================
  // == FFBond class bindings ==================================================
  // ===========================================================================
  py::class_<FFBond>(m, "FFBond")
  .def(py::init<>())
  .def("GetForceConstant", &FFBond::GetForceConstant)
  .def("GetIdealLength", &FFBond::GetIdealLength)
  .def("GetType", &FFBond::GetType)
  .def("GetID", &FFBond::GetID)
  .def("GetLinkedType", &FFBond::GetLinkedType)
  .def("GetForcefield", &FFBond::GetForcefield)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &FFBond::operator bool)
  .def("__str__", &outstream_operator<FFBond>)
  .def("__repr__", &outstream_operator<FFBond>)
  ;
  
  // ===========================================================================
  // == FFAngle class bindings =================================================
  // ===========================================================================
  py::class_<FFAngle>(m, "FFAngle")
  .def(py::init<>())
  .def("GetForceConstant", &FFAngle::GetForceConstant)
  .def("GetIdealAngle", &FFAngle::GetIdealAngle)
  .def("GetType", &FFAngle::GetType)
  .def("GetID", &FFAngle::GetID)
  .def("GetLinkedType", &FFAngle::GetLinkedType)
  .def("GetForcefield", &FFAngle::GetForcefield)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &FFAngle::operator bool)
  .def("__str__", &outstream_operator<FFAngle>)
  .def("__repr__", &outstream_operator<FFAngle>)
  ;
  
  // ===========================================================================
  // == FFDihedral class bindings ==============================================
  // ===========================================================================
  py::class_<FFDihedral>(m, "FFDihedral")
  .def(py::init<>())
  .def("GetPhaseShift", &FFDihedral::GetPhaseShift)
  .def("GetForceConstant", &FFDihedral::GetForceConstant)
  .def("GetMultiplicity", &FFDihedral::GetMultiplicity)
  .def("GetIdealAngle", &FFDihedral::GetIdealAngle)
  .def("GetType", &FFDihedral::GetType)
  .def("GetID", &FFDihedral::GetID)
  .def("GetForcefield", &FFDihedral::GetForcefield)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def(py::self < py::self)
  .def(py::self > py::self)
  .def(py::self <= py::self)
  .def(py::self >= py::self)
  .def("__bool__", &FFDihedral::operator bool)
  .def("__str__", &outstream_operator<FFDihedral>)
  .def("__repr__", &outstream_operator<FFDihedral>)
  ;
  
  // ===========================================================================
  // == Forcefield class bindings ==============================================
  // ===========================================================================
  py::class_<Forcefield>(m, "Forcefield")
  .def(py::init<>())
  .def(py::init<FFFamily, std::string>())
  .def("NewAtomType", &Forcefield::NewAtomType)
  .def("ReserveAtomTypes", &Forcefield::ReserveAtomTypes)
  .def("NewBondType", py::overload_cast<BondType, int32_t, double, double>(&Forcefield::NewBondType))
  .def("LinkBondTypes", &Forcefield::LinkBondTypes)
  .def("ReserveBondTypes", &Forcefield::ReserveBondTypes)
  .def("NewAngleType", py::overload_cast<AngleType, int32_t, double, double>(&Forcefield::NewAngleType))
  .def("LinkAngleTypes", &Forcefield::LinkAngleTypes)
  .def("ReserveAngleTypes", &Forcefield::ReserveAngleTypes)
  .def("NewDihedralType", py::overload_cast<DihedralType, int32_t, double, double, double>(&Forcefield::NewDihedralType))
  .def("NewDihedralType", py::overload_cast<DihedralType, int32_t, double, double>(&Forcefield::NewDihedralType))
  .def("ReserveDihedralTypes", &Forcefield::ReserveDihedralTypes)
  .def("GetAtomType", py::overload_cast<std::string>(&Forcefield::GetAtomType, py::const_))
  .def("GetAtomType", py::overload_cast<int32_t>(&Forcefield::GetAtomType, py::const_))
  .def("NumAtomTypes", &Forcefield::NumAtomTypes)
  .def("GetBondType", py::overload_cast<BondType, int32_t>(&Forcefield::GetBondType, py::const_))
  .def("GetBondType", py::overload_cast<int32_t>(&Forcefield::GetBondType, py::const_))
  .def("NumBondTypes", py::overload_cast<>(&Forcefield::NumBondTypes, py::const_))
  .def("NumBondTypes", py::overload_cast<BondType>(&Forcefield::NumBondTypes, py::const_))
  .def("GetAngleType", py::overload_cast<AngleType, int32_t>(&Forcefield::GetAngleType, py::const_))
  .def("GetAngleType", py::overload_cast<int32_t>(&Forcefield::GetAngleType, py::const_))
  .def("NumAngleTypes", py::overload_cast<>(&Forcefield::NumAngleTypes, py::const_))
  .def("NumAngleTypes", py::overload_cast<AngleType>(&Forcefield::NumAngleTypes, py::const_))
  .def("GetDihedralType", py::overload_cast<DihedralType, int32_t>(&Forcefield::GetDihedralType, py::const_))
  .def("GetDihedralType", py::overload_cast<int32_t>(&Forcefield::GetDihedralType, py::const_))
  .def("NumDihedralTypes", py::overload_cast<>(&Forcefield::NumDihedralTypes, py::const_))
  .def("NumDihedralTypes", py::overload_cast<DihedralType>(&Forcefield::NumDihedralTypes, py::const_))
  .def("GetFamily", &Forcefield::GetFamily)
  .def("GetName", &Forcefield::GetName)
  .def("__bool__", &Forcefield::operator bool)
  .def(py::self == py::self)
  .def(py::self != py::self)
  .def("__str__", &outstream_operator<Forcefield>)
  .def("__repr__", &outstream_operator<Forcefield>)
  ;
  
  // ===========================================================================
  // == Module function bindings ===============================================
  // ===========================================================================
  m.def("GenerateGROMOS54A7", &GenerateGROMOS54A7);
}

