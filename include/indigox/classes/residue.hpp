/*! \file residue.hpp */
#include "../utils/fwd_declares.hpp"

#include <EASTL/vector_set.h>
#include <vector>

#ifndef INDIGOX_CLASSES_RESIDUE_HPP
#define INDIGOX_CLASSES_RESIDUE_HPP

namespace indigox {
  class Residue {
    friend class cereal::access;
    friend class Molecule;

  public:
    using ResidueAtoms = eastl::vector_set<Atom>;

  public:
    enum class Type {
      Unspecified,
      AminoAcid,
      NucleicAcid,
      Sugar,
      Lipid,
      NonSpecific
    };

  private:
    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

  private:
    Residue(const std::vector<Atom> &atoms, const Molecule &mol);

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Residue);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Residue, res);

  public:
    bool HasAtom(const Atom &atom) const;
    /*! \brief Determine the type of the residue.
     *  \details Checks if the residue matches each of the specific residue
     *  types, in order. The first type which matches is the returned type. If
     *  you want to check if a residue is also another type, use the IsType
     *  methods directly.
     */
    Type GetType();
    
    /*! \brief Determine if the residue is an amino acid.
     *  \details The algorithm for determining if a residue is an amino acid is
     *  as follows. First, all potential carboxylic acid and amine atoms are
     *  identified. For carbon atoms, they are classed as potential carboxylic
     *  acid atoms if they are 3 co-ordinate, have a carbonyl bond and another
     *  neighbouring oxygen atom. The additional neighbouring oxygen atom can
     *  also be accounted for by a bond between this residue and the next. For
     *  nitrogen atoms, they are classed as potential amine atoms if they have a
     *  hydrogen bond or a bond to the next residue. This allows for proline to
     *  be identified as an amino acid even though it technically is not.
     *
     *  Once all potential carboxylic and amine atoms are identified, shortest
     *  paths between all pairs of them are checked. If one of these paths
     *  consists of only carbon atoms, the residue is regarded as being an
     *  amino acid.
     */
    bool IsAminoAcid();
    bool IsAlphaAminoAcid();
    bool IsBetaAminoAcid();
    bool IsGammaAminoAcid();
    bool IsDeltaAminoAcid();
    bool IsSugar();
    bool IsLipid();
    bool IsNucleicAcid();
    const ResidueAtoms &GetAtoms() const;

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
  };

  using ResidueType = Residue::Type;
} // namespace indigox

/*
 struct Residue::Impl {
 ResidueType type;
 ResidueAtoms atoms;
 Molecule molecule;
 graph::MolecularGraph residue_graph;

 template <typename Archive>
 void serialise(Archive &archive, const uint32_t version);

 Impl() = default;
 Impl(const std::vector<Atom> &atms, const Molecule &mol);
 };
 */

#endif /* INDIGOX_CLASSES_RESIDUE_HPP */
