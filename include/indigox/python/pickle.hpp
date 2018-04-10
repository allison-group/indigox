#ifndef INDIGOX_PYTHON_PICKLE_HPP
#define INDIGOX_PYTHON_PICKLE_HPP

#include <pybind11/pybind11.h>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>

namespace py = pybind11;

py::tuple PickleAtom(const indigox::Atom a);
indigox::Atom UnpickleAtom(py::tuple& t);

#endif /* INDIGOX_PYTHON_PICKLE_HPP */
