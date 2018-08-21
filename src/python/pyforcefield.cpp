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
  .def("GetForcefield", &IXFFAtom::GetForcefield)
  ;
  
  py::class_<IXFFBond, FFBond>(m, "FFBond")
  // No constructor
  .def("GetForceConstant", &IXFFBond::GetForceConstant)
  .def("GetIdealLength", &IXFFBond::GetIdealLength)
  .def("GetType", &IXFFBond::GetType)
  .def("GetID", &IXFFBond::GetID)
  .def("GetLinkedType", &IXFFBond::GetLinkedType)
  .def("GetForcefield", &IXFFBond::GetForcefield)
  ;
  
  py::class_<IXFFAngle, FFAngle>(m, "FFAngle")
  // No constructor
  .def("GetForceConstant", &IXFFAngle::GetForceConstant)
  .def("GetIdealAngle", &IXFFAngle::GetIdealAngle)
  .def("GetType", &IXFFAngle::GetType)
  .def("GetID", &IXFFAngle::GetID)
  .def("GetLinkedType", &IXFFAngle::GetLinkedType)
  .def("GetForcefield", &IXFFAngle::GetForcefield)
  ;
  
  py::class_<IXFFDihedral, FFDihedral>(m, "FFDihedral")
  // No constructor
  .def("GetPhaseShift", &IXFFDihedral::GetPhaseShift)
  .def("GetForceConstant", &IXFFDihedral::GetForceConstant)
  .def("GetMultiplicity", &IXFFDihedral::GetMultiplicity)
  .def("GetIdealAngle", &IXFFDihedral::GetIdealAngle)
  .def("GetType", &IXFFDihedral::GetType)
  .def("GetID", &IXFFDihedral::GetID)
  .def("GetForcefield", &IXFFDihedral::GetForcefield)
  ;
  
  py::class_<IXForcefield, Forcefield>(m, "Forcefield")
  .def(py::init<FFFamily, string_>())
  .def("GetFamily", &IXForcefield::GetFamily)
  .def("GetName", &IXForcefield::GetName)
  // Atoms
  .def("NewAtomType", &IXForcefield::NewAtomType)
  .def("GetAtomType", py::overload_cast<string_>(&IXForcefield::GetAtomType, py::const_))
  .def("GetAtomType", py::overload_cast<int_>(&IXForcefield::GetAtomType, py::const_))
  .def("ReserveAtomTypes", &IXForcefield::ReserveAtomTypes)
  .def("NumAtomTypes", &IXForcefield::NumAtomTypes)
  // Bonds
  .def("NewBondType", py::overload_cast<BondType, int_, float_, float_>(&IXForcefield::NewBondType))
  .def("GetBondType", py::overload_cast<BondType, int_>(&IXForcefield::GetBondType, py::const_))
  .def("GetBondType", py::overload_cast<int_>(&IXForcefield::GetBondType, py::const_))
  .def("LinkBondTypes", &IXForcefield::LinkBondTypes)
  .def("ReserveBondTypes", &IXForcefield::ReserveBondTypes)
  .def("NumBondTypes", py::overload_cast<>(&IXForcefield::NumBondTypes, py::const_))
  .def("NumBondTypes", py::overload_cast<BondType>(&IXForcefield::NumBondTypes, py::const_))
  // Angles
  .def("NewAngleType", py::overload_cast<AngleType, int_, float_, float_>(&IXForcefield::NewAngleType))
  .def("GetAngleType", py::overload_cast<AngleType, int_>(&IXForcefield::GetAngleType, py::const_))
  .def("GetAngleType", py::overload_cast<int_>(&IXForcefield::GetAngleType, py::const_))
  .def("LinkAngleTypes", &IXForcefield::LinkAngleTypes)
  .def("ReserveAngleTypes", &IXForcefield::ReserveAngleTypes)
  .def("NumAngleTypes", py::overload_cast<>(&IXForcefield::NumAngleTypes, py::const_))
  .def("NumAngleTypes", py::overload_cast<AngleType>(&IXForcefield::NumAngleTypes, py::const_))
  // Dihedrals
  .def("NewDihedralType", py::overload_cast<DihedralType, int_, float_, float_, float_>(&IXForcefield::NewDihedralType))
  .def("NewDihedralType", py::overload_cast<DihedralType, int_, float_, float_>(&IXForcefield::NewDihedralType))
  .def("GetDihedralType", py::overload_cast<DihedralType, int_>(&IXForcefield::GetDihedralType, py::const_))
  .def("GetDihedralType", py::overload_cast<int_>(&IXForcefield::GetDihedralType, py::const_))
  .def("ReserveDihedralTypes", &IXForcefield::ReserveDihedralTypes)
  .def("NumDihedralTypes", py::overload_cast<>(&IXForcefield::NumDihedralTypes, py::const_))
  .def("NumDihedralTypes", py::overload_cast<DihedralType>(&IXForcefield::NumDihedralTypes, py::const_))
  ;
}

