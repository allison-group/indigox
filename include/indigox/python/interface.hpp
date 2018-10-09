#include <pybind11/pybind11.h>

#ifndef INDIGOX_PYTHON_INTERFACE_HPP
#define INDIGOX_PYTHON_INTERFACE_HPP

#include "../utils/fwd_declares.hpp"

void GeneratePyAngle(pybind11::module& m);
void GeneratePyAtom(pybind11::module& m);
void GeneratePyBond(pybind11::module& m);
void GeneratePyDihedral(pybind11::module& m);
void GeneratePyMolecule(pybind11::module& m);
void GeneratePyMolecularGraph(pybind11::module& m);
void GeneratePyCondensedMolecularGraph(pybind11::module& m);
void GeneratePyPeriodicTable(pybind11::module& m);
void GeneratePyElement(pybind11::module& m);
void GeneratePyElectronAssignmentGraph(pybind11::module& m);
void GeneratePyElectronAssigner(pybind11::module& m);
void GeneratePyForcefield(pybind11::module& m);
void GeneratePyAthenaeum(pybind11::module& m);
void GenerateOpaqueContainers(pybind11::module& m);

#endif /* INDIGOX_PYTHON_INTERFACE_HPP */
