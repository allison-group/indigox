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
  
  py::enum_<ForcefieldFamily>(m, "ForcefieldFamily")
  .value("Empty", ForcefieldFamily::Empty)
  .value("GROMOS", ForcefieldFamily::GROMOS)
  .value("CHARMM", ForcefieldFamily::CHARMM)
  .value("AMBER", ForcefieldFamily::AMBER)
  .value("Other", ForcefieldFamily::Other)
  ;
  
  py::enum_<BondFunctionType>(m, "BondFunctionType")
  .value("Empty", BondFunctionType::Empty)
  .value("Harmonic", BondFunctionType::Harmonic)
  .value("Quartic", BondFunctionType::Quartic)
  .value("Morse", BondFunctionType::Morse)
  .value("Cubic", BondFunctionType::Cubic)
  .value("FENE", BondFunctionType::FENE)
  ;
  
  py::enum_<AngleFunctionType>(m, "AngleFunctionType")
  .value("Empty", AngleFunctionType::Empty)
  .value("Harmonic", AngleFunctionType::Harmonic)
  .value("CosineHarmonic", AngleFunctionType::CosineHarmonic)
  .value("UreyBradley", AngleFunctionType::UreyBradley)
  .value("Quartic", AngleFunctionType::Quartic)
  ;
  
  py::enum_<DihedralFunctionType>(m, "DihedralFunctionType")
  .value("Empty", DihedralFunctionType::Empty)
  .value("Proper", DihedralFunctionType::Proper)
  .value("Improper", DihedralFunctionType::Improper)
  .value("RyckaertBellemans", DihedralFunctionType::RyckaertBellemans)
  .value("PeriodicImproper", DihedralFunctionType::PeriodicImproper)
  .value("Fourier", DihedralFunctionType::Fourier)
  .value("Restricted", DihedralFunctionType::Restricted)
  ;
  
  py::class_<IXFFAtomType, FFAtomType>(m, "FFAtomType")
  // No constructor
  .def("GetID", &IXFFAtomType::GetID)
  .def("GetName", &IXFFAtomType::GetName)
  ;
  
  py::class_<IXFFBondType, FFBondType>(m, "FFBondType")
  // No constructor
  .def("GetForceConstant", &IXFFBondType::GetForceConstant)
  .def("GetIdealLength", &IXFFBondType::GetIdealLength)
  .def("GetType", &IXFFBondType::GetType)
  .def("GetID", &IXFFBondType::GetID)
  ;
  
  py::class_<IXFFAngleType, FFAngleType>(m, "FFAngleType")
  // No constructor
  .def("GetForceConstant", &IXFFAngleType::GetForceConstant)
  .def("GetIdealAngle", &IXFFAngleType::GetIdealAngle)
  .def("GetType", &IXFFAngleType::GetType)
  .def("GetID", &IXFFAngleType::GetID)
  ;
  
  py::class_<IXFFDihedralType, FFDihedralType>(m, "FFDihedralType")
  // No constructor
  .def("GetPhaseShift", &IXFFDihedralType::GetPhaseShift)
  .def("GetForceConstant", &IXFFDihedralType::GetForceConstant)
  .def("GetMultiplicity", &IXFFDihedralType::GetMultiplicity)
  .def("GetIdealAngle", &IXFFDihedralType::GetIdealAngle)
  .def("GetType", &IXFFDihedralType::GetType)
  .def("GetID", &IXFFDihedralType::GetID)
  ;
  
  py::class_<IXForcefield, Forcefield>(m, "Forcefield")
  .def(py::init<ForcefieldFamily, string_>())
  .def("ToString", &IXForcefield::ToString)
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

