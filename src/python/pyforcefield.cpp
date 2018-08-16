#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <indigox/python/interface.hpp>
#include <indigox/python/pickle.hpp>

#include <indigox/classes/forcefield.hpp>

namespace py = pybind11;

void GeneratePyForcefield(py::module& m) {
  using namespace indigox;
  
  m.def("GenerateGROMOS54A7", &GenerateGROMOS54A7);
  
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
  
  py::class_<IXFFAtom, FFAtom>(m, "FFAtom")
  // No constructor
  .def("GetID", &IXFFAtom::GetID)
  .def("GetName", &IXFFAtom::GetName)
  ;
  
  py::class_<IXFFBond, FFBond>(m, "FFBond")
  // No constructor
  .def("GetForceConstant", &IXFFBond::GetForceConstant)
  .def("GetIdealLength", &IXFFBond::GetIdealLength)
  .def("GetType", &IXFFBond::GetType)
  .def("GetID", &IXFFBond::GetID)
  ;
  
  py::class_<IXFFAngle, FFAngle>(m, "FFAngle")
  // No constructor
  .def("GetForceConstant", &IXFFAngle::GetForceConstant)
  .def("GetIdealAngle", &IXFFAngle::GetIdealAngle)
  .def("GetType", &IXFFAngle::GetType)
  .def("GetID", &IXFFAngle::GetID)
  ;
  
  py::class_<IXFFDihedral, FFDihedral>(m, "FFDihedral")
  // No constructor
  .def("GetPhaseShift", &IXFFDihedral::GetPhaseShift)
  .def("GetForceConstant", &IXFFDihedral::GetForceConstant)
  .def("GetMultiplicity", &IXFFDihedral::GetMultiplicity)
  .def("GetIdealAngle", &IXFFDihedral::GetIdealAngle)
  .def("GetType", &IXFFDihedral::GetType)
  .def("GetID", &IXFFDihedral::GetID)
  ;
  
  py::class_<IXForcefield, Forcefield>(m, "Forcefield")
  .def(py::init<FFFamily, string_>())
  .def("ReserveAtomTypes", &IXForcefield::ReserveAtomTypes)
  .def("ReserveBondTypes", &IXForcefield::ReserveBondTypes)
  .def("ReserveAngleTypes", &IXForcefield::ReserveAngleTypes)
  .def("ReserveDihedralTypes", &IXForcefield::ReserveDihedralTypes)
  .def("GetFamily", &IXForcefield::GetFamily)
  .def("GetName", &IXForcefield::GetName)
  .def("NewAtomType", &IXForcefield::NewAtomType)
  .def("GetAtomType", py::overload_cast<string_>(&IXForcefield::GetAtomType, py::const_))
  .def("GetAtomType", py::overload_cast<int_>(&IXForcefield::GetAtomType, py::const_))
  .def("NewHarmonicBondType", &IXForcefield::NewHarmonicBondType)
  .def("NewQuarticBondType", &IXForcefield::NewQuarticBondType)
  .def("GetBondType", &IXForcefield::GetBondType)
  .def("GetHarmonicBondType", &IXForcefield::GetHarmonicBondType)
  .def("GetQuarticBondType", &IXForcefield::GetQuarticBondType)
  .def("NewHarmonicAngleType", &IXForcefield::NewHarmonicAngleType)
  .def("NewCosineHarmonicAngleType", &IXForcefield::NewCosineHarmonicAngleType)
  .def("GetAngleType", &IXForcefield::GetAngleType)
  .def("GetHarmonicAngletype", &IXForcefield::GetHarmonicAngleType)
  .def("GetCosineHarmonicAngleType", &IXForcefield::GetCosineHarmonicAngleType)
  .def("NewProperDihedralType", &IXForcefield::NewProperDihedralType)
  .def("NewImproperDihedralType", &IXForcefield::NewImproperDihedralType)
  .def("GetDihedralType", &IXForcefield::GetDihedralType)
  .def("GetImproperDihedralType", &IXForcefield::GetImproperDihedralType)
  .def("GetProperDihedralType", &IXForcefield::GetProperDihedralType)
  .def("NumAtomTypes", &IXForcefield::NumAtomTypes)
  .def("NumBondTypes", &IXForcefield::NumBondTypes)
  .def("NumAngleTypes", &IXForcefield::NumAngleTypes)
  .def("NumDihedralTypes", &IXForcefield::NumDihedralTypes)
  ;
}

