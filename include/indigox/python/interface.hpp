#include <sstream>

#include <pybind11/pybind11.h>


#ifndef INDIGOX_PYTHON_INTERFACE_HPP
#define INDIGOX_PYTHON_INTERFACE_HPP

#include "../utils/fwd_declares.hpp"

void GeneratePyMolecule(pybind11::module& m);
void GeneratePyPeriodicTable(pybind11::module& m);
void GeneratePyGraphs(pybind11::module& m);
void GeneratePyForcefield(pybind11::module& m);
void GeneratePyAthenaeum(pybind11::module& m);
void GeneratePyGraphAlgorithms(pybind11::module& m);

void GeneratePyElectronAssigner(pybind11::module& m);

void GenerateOpaqueContainers(pybind11::module& m);

template<typename T>
std::string outstream_operator(const T& t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

#endif /* INDIGOX_PYTHON_INTERFACE_HPP */
