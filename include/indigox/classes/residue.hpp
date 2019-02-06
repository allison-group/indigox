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
    Type GetType();
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
