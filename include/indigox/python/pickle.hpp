#ifndef INDIGOX_PYTHON_PICKLE_HPP
#define INDIGOX_PYTHON_PICKLE_HPP

#include <pybind11/pybind11.h>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>

namespace py = pybind11;

py::tuple PickleAtom(const indigox::Atom atom);
indigox::Atom UnpickleAtom(py::tuple& t);

py::tuple PickleBond(const indigox::Bond bond);
indigox::Bond UnpickleBond(py::tuple& t);

py::tuple PicklePeriodicTable(const indigox::PeriodicTable table);
indigox::PeriodicTable UnpicklePeriodicTable(py::tuple& t);

#endif /* INDIGOX_PYTHON_PICKLE_HPP */
