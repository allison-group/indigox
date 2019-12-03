#include "../utils/fwd_declares.hpp"
#include "../utils/quad.hpp"
#include "../utils/triple.hpp"
#include <boost/dynamic_bitset_fwd.hpp>

#include <algorithm>
#include <map>
#include <vector>

#ifndef INDIGOX_CLASSES_ATHENAEUM_HPP
#define INDIGOX_CLASSES_ATHENAEUM_HPP

namespace indigox {

  /*! \brief Fragment class for CherryPicker parameterisation algorithm.
   */
  class Fragment {
    friend class cereal::access;
    friend class Athenaeum;

  public:
    /*! \brief Type of overlapping vertex.
     *  \details Each vertex of an overlap is given a type. This is to enable
     *  more detailed isomorphism matching methods to be used when matching. */
    enum class OverlapType { GenericOverlap };

    //! \brief Type used to match an overlap vertex to its type of overlap.
    using OverlapVertex = std::pair<OverlapType, graph::CMGVertex>;
    using Vert = graph::MGVertex;
    using AtmType = graph::MGVertex;
    using BndType = std::pair<AtmType, AtmType>;
    using AngType = stdx::triple<AtmType>;
    using DhdType = stdx::quad<AtmType>;

  private:
    template <class Archive>
    void serialise(Archive &archive, const uint32_t version);

  public:
    //! \cond ignored
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Fragment);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Fragment, frag);
    //! \endcond

    /*! \brief Normal constructor
     *  \details Constructs a fragment of the given graph from the combined
     *  fragment and overlap vertices.
     *  \param G the graph to generate fragment from.
     *  \param frag the vertices of G that are part of the fragment.
     *  \param overlap the vertices of G that are part of the overlap region. */
    Fragment(const graph::MolecularGraph &G, std::vector<graph::MGVertex> &frag,
             std::vector<graph::MGVertex> &overlap);

    Fragment(const Molecule &mol, std::vector<Atom> &frag,
             std::vector<Atom> &overlap);

    const graph::CondensedMolecularGraph &GetGraph() const;
    const std::vector<graph::CMGVertex> &GetFragment() const;
    const boost::dynamic_bitset<> &GetSupersets() const;
    size_t Size() const;
    const std::vector<OverlapVertex> &GetOverlap() const;
    bool IsFragmentVertex(const graph::CMGVertex &v) const;
    bool IsOverlapVertex(const graph::CMGVertex &v) const;
    const std::vector<AtmType> &GetAtoms() const;
    const std::vector<BndType> &GetBonds() const;
    const std::vector<AngType> &GetAngles() const;
    const std::vector<DhdType> &GetDihedrals() const;

  private:
    struct FragmentData;
    std::shared_ptr<FragmentData> m_data;
  };

  /*! \brief Athenaeum class for fragment storage in CherryPicker algorithm.
   */
  class Athenaeum {
    friend class cereal::access;

  public:
    enum class Settings : uint8_t {
      FragmentCycles,
      SelfConsistent,
      BoolCount,
      MoleculeSizeLimit,
      MinimumFragmentSize,
      MaximumFragmentSize,
      OverlapLength,
      CycleSize,
      IntCount
    };

    // Settings
    bool GetBool(Settings param);
    void SetBool(Settings param);
    void UnsetBool(Settings param);

    int32_t GetInt(Settings param);
    void SetInt(Settings param, int32_t value);

    void DefaultSettings();
    // End Settings

    using FragContain = std::vector<Fragment>;
    using MoleculeFragments = std::map<Molecule, FragContain>; //map of molecule to vectors of fragments

  private:
    template <class Archive>
    void serialise(Archive &archive, const uint32_t version);

  public:
    //! \cond ignored
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Athenaeum);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Athenaeum, ath);
    //! \endcond

    Athenaeum(const Forcefield &ff);
    Athenaeum(const Forcefield &ff, int32_t overlap);
    Athenaeum(const Forcefield &ff, int32_t overlap, int32_t ring_overlap);

    size_t NumFragments() const;
    size_t NumFragments(const Molecule &mol) const;

    const MoleculeFragments &GetFragments() const;
    const FragContain &GetFragments(const Molecule &mol) const;
    bool HasFragments(const Molecule &mol) const;

    const Forcefield &GetForcefield() const;
    //    bool CheckSelfConsistent();

    /*! \brief Adds the given fragment.
     *  \returns if the fragment was added or not. */
    bool AddFragment(const Fragment &frag);

    /*! \brief Determines all the fragments of a molecule and adds them.
     *  \returns the number of fragments added. */
    size_t AddAllFragments(const Molecule &mol);

  private:
    void SortAndMask(const Molecule &mol);

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
  };

  void SaveAthenaeum(const Athenaeum &ath, const std::string& path);
  Athenaeum LoadAthenaeum(std::string path);

} // namespace indigox

#endif /* INDIGOX_CLASSES_ATHENAEUM_HPP */
