#include "../utils/fwd_declares.hpp"

#include <EASTL/vector_map.h>
#include <cstdint>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>

#ifndef INDIGOX_PYTHON_PYOPAQUECONTAINERS_HPP
#define INDIGOX_PYTHON_PYOPAQUECONTAINERS_HPP

// PYBIND11_MAKE_OPAQUE(eastl::vector_set<indigox::Element>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::CMGVertex>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::CMGEdge>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::MGVertex>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::MGEdge>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<indigox::graph::MGEdge>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<indigox::graph::CMGEdge>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<indigox::graph::MGVertex>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<indigox::graph::CMGVertex>>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::CondensedMolecularGraph>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::graph::MolecularGraph>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::Atom>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::Bond>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::Angle>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::Dihedral>);

PYBIND11_MAKE_OPAQUE(std::vector<indigox::ParamAtom>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::ParamBond>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::ParamAngle>);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::ParamDihedral>);
using mapffatomsize = eastl::vector_map<indigox::FFAtom, size_t>;
using mapffbondsize = eastl::vector_map<indigox::FFBond, size_t>;
using mapffanglesize = eastl::vector_map<indigox::FFAngle, size_t>;
using mapffdihedralsize =
    eastl::vector_map<std::vector<indigox::FFDihedral>, size_t>;
PYBIND11_MAKE_OPAQUE(mapffatomsize);
PYBIND11_MAKE_OPAQUE(mapffbondsize);
PYBIND11_MAKE_OPAQUE(mapffanglesize);
PYBIND11_MAKE_OPAQUE(std::vector<indigox::FFDihedral>);
PYBIND11_MAKE_OPAQUE(mapffdihedralsize);

PYBIND11_MAKE_OPAQUE(std::vector<indigox::Fragment>);
using mapmolvecfrag =
    std::map<indigox::Molecule, std::vector<indigox::Fragment>>;
PYBIND11_MAKE_OPAQUE(mapmolvecfrag);

#endif /* INDIGOX_PYTHON_PYOPAQUECONTAINERS_HPP */
