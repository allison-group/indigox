#ifndef INDIGOX_CLASSES_MOLECULE_IMPL_HPP
#define INDIGOX_CLASSES_MOLECULE_IMPL_HPP

#include "../graph/condensed.hpp"
#include "../graph/molecular.hpp"
#include "../utils/fwd_declares.hpp"
#include "angle.hpp"
#include "atom.hpp"
#include "bond.hpp"
#include "dihedral.hpp"
#include "forcefield.hpp"
#include "molecule.hpp"
#include "periodictable.hpp"
#include "residue.hpp"

#include <bitset>

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
    int32_t charge_group_id;
    int32_t residue_id;
    int32_t implicit_hydrogens;
    Eigen::Vector3d position;
    std::string name;
    std::string residue_name;
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
  // == RESIDUE IMPLEMENTATION =============================================
  // =======================================================================

  struct Residue::Impl {
    ResidueType type;
    ResidueAtoms atoms;
    Molecule molecule;
    graph::MolecularGraph residue_graph;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

    Impl() = default;
    Impl(const std::vector<Atom> &atms, const Molecule &mol);

    bool AminoAcidTest();

    void DetermineType();
  };

  // =======================================================================
  // == MOLECULE IMPLEMENTATION ============================================
  // =======================================================================

  enum class CalculatedData : uint8_t {
    Formula,
    AnglePerception,
    DihedralPerception,
    ResiduePerception,
    ChargeGroupDetermination,
    CondensedGraph,
    Number
  };

  struct Molecule::Impl {
    std::string name;
    int64_t next_unique_id;
    int32_t molecular_charge;
    MoleculeAtoms atoms;
    MoleculeBonds bonds;
    MoleculeAngles angles;
    MoleculeDihedrals dihedrals;
    MoleculeResidues residues;
    Forcefield forcefield;
    graph::MolecularGraph molecular_graph;
    graph::CondensedMolecularGraph condensed_molecular_graph;

    std::bitset<static_cast<uint8_t>(CalculatedData::Number)> calculated_data;

    // Cached variables
    std::string cached_formula;

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t);

    Impl() = default;
    Impl(std::string n);

    // Return the index in a's bonds if found, -1 if not found
    int64_t FindBond(const Atom &a, const Atom &b) const;
    // Return the index in b's angles if found, -1 if not found
    int64_t FindAngle(const Atom &a, const Atom &b, const Atom &c) const;
    // Return the index in a's dihedrals if found, -1 if not found
    int64_t FindDihedral(const Atom &a, const Atom &b, const Atom &c,
                         const Atom &d) const;

    inline bool Test(CalculatedData dat) const {
      return calculated_data.test(static_cast<uint8_t>(dat));
    }
    inline void Set(CalculatedData dat) {
      calculated_data.set(static_cast<uint8_t>(dat));
    }
    inline void Reset(CalculatedData dat) {
      calculated_data.reset(static_cast<uint8_t>(dat));
    }
    inline void Reset() { calculated_data.reset(); }
  };
} // namespace indigox

#endif /* molecule_impl_hpp */
