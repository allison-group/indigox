#include <pybind11/pybind11.h>

#ifndef INDIGOX_PYTHON_INTERFACE_HPP
#define INDIGOX_PYTHON_INTERFACE_HPP

/// @todo add opaque set<string> stuff so can add to string
// PYBIND11_MAKE_OPAQUE(std::set<indigox::String>)


void GenerateOptions(pybind11::module& m);
void GeneratePyAngle(pybind11::module& m);
void GeneratePyAtom(pybind11::module& m);
void GeneratePyBond(pybind11::module& m);
void GeneratePyDihedral(pybind11::module& m);
void GeneratePyMolecule(pybind11::module& m);
void GeneratePyMolecularGraph(pybind11::module& m);
void GeneratePyPeriodicTable(pybind11::module& m);
void GeneratePyElement(pybind11::module& m);
void GeneratePyElectronAssignmentGraph(pybind11::module& m);



#endif /* INDIGOX_PYTHON_INTERFACE_HPP */
