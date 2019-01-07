#include <indigox/algorithm/electron_assignment/assigner.hpp>
#include <indigox/algorithm/electron_assignment/astar_optimisation.hpp>
#include <indigox/algorithm/electron_assignment/local_optimisation.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/python/interface.hpp>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void GeneratePyElectronAssigner(pybind11::module &m) {
  using namespace indigox::algorithm;
  using Settings = IXElectronAssigner::Settings;
  using LOSetters = IXLocalOptimisation::Settings;
  using ASSettings = IXAStarOptimisation::Settings;
  using Algorithm = IXElectronAssigner::AssignerAlgorithm;

  //  ElectronAssigner class
  py::class_<IXElectronAssigner, ElectronAssigner> PyAssigner(
      m, "ElectronAssigner");
  PyAssigner.def("Run", &IXElectronAssigner::Run)
      .def("GetAssignmentGraph", &IXElectronAssigner::GetAssignmentGraph)
      .def("LoadScoreTable", &IXElectronAssigner::LoadScoreTable)
      .def("ApplyAssignment", &IXElectronAssigner::ApplyAssignment)
      .def("GetOptimalCount", &IXElectronAssigner::GetOptimalCount)
      .def("GetOptimisedScore", &IXElectronAssigner::GetOptimisedScore);
  m.def("CreateElectronAssigner", &CreateElectronAssigner);

  //  Algorithm enum
  py::enum_<Algorithm>(m, "AssignerAlgorithm")
      .value("LocalOptimisation", Algorithm::LocalOptimisation)
      .value("AStar", Algorithm::AStar)
      .value("FPT", Algorithm::FPT);

  //  ElectronAssigner Settings
  py::class_<Settings, std::shared_ptr<Settings>> PySettings(PyAssigner,
                                                             "Settings");
  PySettings.def_readwrite_static("Algorithm", &Settings::Algorithm)
      .def_readwrite_static("ElectronPairs", &Settings::ElectronPairs)
      .def_readwrite_static("ChargedCarbon", &Settings::ChargedCarbon)
      .def_readwrite_static("Preassign", &Settings::Preassign)
      .def_readwrite_static("ScoreFile", &Settings::ScoreFile)
      .def_readwrite_static("Infinity", &Settings::Infinity)
      .def_readwrite_static("MaxBondOrder", &Settings::MaxBondOrder)
      .def_readwrite_static("MaxChargeMagnitude", &Settings::MaxChargeMagnitude)
      .def_readwrite_static("MaxNumResults", &Settings::MaxNumResults)
      .def_readwrite_static("AllowedElements", &Settings::AllowedElements);

  //  LocalOptimisation Settings
  py::class_<LOSetters, std::shared_ptr<LOSetters>>(PySettings,
                                                    "LocalOptimisation")
      .def_readwrite_static("OptimistationLevel", &LOSetters::OPTIMISE_LEVEL)
      .def_readwrite_static("Timeout", &LOSetters::TIMEOUT)
      .def_readwrite_static("UseCache", &LOSetters::USE_CACHE);

  //  AStar settings
  py::class_<ASSettings, std::shared_ptr<ASSettings>>(PySettings, "AStar")
      .def_readwrite_static("MemoryLimit", &ASSettings::MEMORY_LIMIT)
      .def_readwrite_static("Heuristic", &ASSettings::HEURISTIC);
  //  AStar Heuristic enum
  py::enum_<IXAStarOptimisation::Heuristic>(PySettings, "AStarHeuristic")
      .value("Promiscuous", IXAStarOptimisation::Heuristic::Promiscuous)
      .value("Abstemious", IXAStarOptimisation::Heuristic::Abstemious);
}
