// This file contains forward declarations of all classes used in indigox as
// well as their shared_ptr and weak_ptr counterparts, as needed.
#ifndef INDIGOX_UTILS_FWD_DECLARES_HPP
#define INDIGOX_UTILS_FWD_DECLARES_HPP
#include <cstdint>
#include <memory>

// Serialisation related stuff, using the cereal library
namespace cereal {
  class access;
  class PortableBinaryInputArchive;
  class PortableBinaryOutputArchive;
  class JSONInputArchive;
  class JSONOutputArchive;
  template <class T> class construct;
} // namespace cereal

#define INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(class_name)                       \
  class_name() = default;                                                      \
  class_name(const class_name &) = default;                                    \
  class_name(class_name &&) = default;                                         \
  class_name &operator=(const class_name &) = default;                         \
  class_name &operator=(class_name &&) = default;                              \
  ~class_name() = default

#define INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(class_name, short_name)          \
  bool operator==(const class_name &short_name) const;                         \
  bool operator<(const class_name &short_name) const;                          \
  bool operator>(const class_name &short_name) const;                          \
  inline bool operator!=(const class_name &short_name) const {                 \
    return !(operator==(short_name));                                          \
  }                                                                            \
  inline bool operator<=(const class_name &short_name) const {                 \
    return !(operator>(short_name));                                           \
  }                                                                            \
  inline bool operator>=(const class_name &short_name) const {                 \
    return !(operator<(short_name));                                           \
  }                                                                            \
  inline operator bool() const {                                               \
    return bool(m_data);                                                       \
  }                                                                            \
  friend std::ostream &operator<<(std::ostream &os,                            \
                                  const class_name &short_name)

namespace indigox {

  using State = uint32_t;

  // Molecule related
  class Molecule;
  class Atom;
  class Bond;
  class Angle;
  class Dihedral;
  class Element;
  class PeriodicTable;

  // CherryPicker Related
  class ParamMolecule;
  class ParamAtom;
  class ParamBond;
  class ParamAngle;
  class ParamDihedral;

  // Forcefield related
  class Forcefield;
  class FFAtom;
  class FFBond;
  class FFAngle;
  class FFDihedral;

  // Athenaeum related
  class Fragment;
  using sFragment = std::shared_ptr<Fragment>;
  using wFragment = std::weak_ptr<Fragment>;

  class Athenaeum;
  using sAthenaeum = std::shared_ptr<Athenaeum>;
  using wAthenaeum = std::weak_ptr<Athenaeum>;

  namespace algorithm {
    struct access;

    class CherryPicker;

    class IXElectronAssigner;
    using ElectronAssigner = std::shared_ptr<IXElectronAssigner>;
    using _ElectronAssigner = std::weak_ptr<IXElectronAssigner>;
  } // namespace algorithm

  namespace graph {
    // Base Graph
    template <class V, class E, class S, class D, class VP, class EP>
    class BaseGraph;
    struct Directed;
    struct Undirected;
    struct GraphLabel;

    // AssignmentGraph
    class IXAssignmentGraph;
    using AssignmentGraph = std::shared_ptr<IXAssignmentGraph>;
    using _AssignmentGraph = std::weak_ptr<IXAssignmentGraph>;

    class IXAGVertex;
    using AGVertex = std::shared_ptr<IXAGVertex>;
    using _AGVertex = std::weak_ptr<IXAGVertex>;

    // MolecularGraph
    class MolecularGraph;
    class MGVertex;
    class MGEdge;

    // CondensedMolecularGraph
    class CondensedMolecularGraph;
    class CMGVertex;
    class CMGEdge;
  } // namespace graph

  namespace test {
    struct TestForcefield;
    struct TestFFAtom;
    struct TestFFBond;
    struct TestFFAngle;
    struct TestFFDihedral;

    struct TestMolecule;
    struct TestAtom;
    struct TestBond;
    struct TestAngle;
    struct TestDihedral;
    struct TestElement;
    struct TestPeriodicTable;

    struct TestAssignmentGraph;
    struct TestCondensedMolecularGraph;
    struct TestCondensedVertex;
    struct TestCondensedEdge;
    struct TestMolecularGraph;
    struct TestMolecularVertex;
    struct TestMolecularEdge;
  } // namespace test
} // namespace indigox

#endif
