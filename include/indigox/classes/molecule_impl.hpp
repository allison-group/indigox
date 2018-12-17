#ifndef INDIGOX_CLASSES_MOLECULE_IMPL_HPP
#define INDIGOX_CLASSES_MOLECULE_IMPL_HPP

#include "../utils/fwd_declares.hpp"
#include "angle.hpp"
#include "atom.hpp"
#include "bond.hpp"
#include "dihedral.hpp"
#include "molecule.hpp"
#include "../graph/molecular.hpp"

namespace indigox {

  // =======================================================================
  // == ANGLE IMPLEMENTATION ===============================================
  // =======================================================================

  struct Angle::Impl {
    AngleAtoms atoms;
    Molecule molecule;
    int64_t tag;
    int64_t unique_id;
    FFAngle forcefield_type;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t);

    Impl() = default;
    Impl(const Atom &a, const Atom &b, const Atom &c, const Molecule &mol);
  };

  // =======================================================================
  // == ATOM IMPLEMENTATION ================================================
  // =======================================================================

  struct Atom::Impl {
    Molecule molecule;
    Element element;
    int32_t formal_charge;
    int64_t tag;
    int64_t unique_id;
    int32_t implicit_hydrogens;
    Eigen::Vector3d position;
    std::string name;
    double partial_charge;
    Stereo stereochemistry;
    FFAtom forcefield_type;
    Atom::AtomBonds bonds;
    Atom::AtomAngles angles;
    Atom::AtomDihedrals dihedrals;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t);

    Impl() = default;
    Impl(const Molecule &mol, const Element &elem, double x, double y, double z,
         std::string n);
  };

  // =======================================================================
  // == BOND IMPLEMENTATION ================================================
  // =======================================================================

  struct Bond::Impl {
    BondAtoms atoms;
    Molecule molecule;
    int64_t tag;
    int64_t unique_id;
    BondOrder order;
    BondStereo stereochemistry;
    FFBond forcefield_type;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

    Impl() = default;
    Impl(const Atom &a, const Atom &b, const Molecule &mol, BondOrder o);
  };

  // =======================================================================
  // == DIHEDRAL IMPLEMENTATION ============================================
  // =======================================================================

  struct Dihedral::Impl {
    DihedralAtoms atoms;
    Molecule molecule;
    int64_t tag;
    int64_t unique_id;
    DihedralTypes forcefield_types;
    int32_t priority;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

    Impl() = default;
    Impl(const Atom &a, const Atom &b, const Atom &c, const Atom &d,
         const Molecule &molecule);
  };

  // =======================================================================
  // == MOLECULE IMPLEMENTATION ============================================
  // =======================================================================

  struct Molecule::Impl {
    std::string name;
    int64_t next_unique_id;
    int32_t molecular_charge;
    MoleculeAtoms atoms;
    MoleculeBonds bonds;
    MoleculeAngles angles;
    MoleculeDihedrals dihedrals;
    Forcefield forcefield;
    graph::MolecularGraph molecular_graph;
    State modification_state;
    bool frozen;

    // Cached variables
    std::string cached_formula;
    State cached_formula_state;
    State angle_percieved_state;
    State dihedral_percieved_state;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t);

    Impl() = default;
    Impl(std::string n, const Molecule &mol);
  };
} // namespace indigox

#endif /* molecule_impl_hpp */
