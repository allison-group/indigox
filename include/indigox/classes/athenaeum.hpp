#include <algorithm>
#include <map>
#include <vector>

#include "../utils/fwd_declares.hpp"
#include "../utils/triple.hpp"
#include "../utils/quad.hpp"
#include "forcefield.hpp"
#include "molecule.hpp"
#include "../graph/condensed.hpp"

#ifndef INDIGOX_CLASSES_ATHENAEUM_HPP
#define INDIGOX_CLASSES_ATHENAEUM_HPP

namespace indigox {
  
  class Fragment {
    friend class cereal::access;
  public:
    /*! \brief Type of overlapping vertex.
     *  \details Each vertex of an overlap is given a type. This is to enable
     *  more detailed isomorphism matching methods to be used when matching. */
    enum class OverlapType {
      GenericOverlap
    };
    
    //! \brief Type used to match an overlap vertex to its type of overlap.
    using OverlapVertex = std::pair<OverlapType, graph::CMGVertex>;
    using Vert = graph::MGVertex;
    using AtmType = graph::MGVertex;
    using BndType = std::pair<AtmType, AtmType>;
    using AngType = stdx::triple<AtmType>;
    using DhdType = stdx::quad<AtmType>;
    
  private:
    template <class Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    Fragment();
    
    /*! \brief Normal constructor
     *  \details Constructs a fragment of the given graph from the combined
     *  fragment and overlap vertices.
     *  \param G the graph to generate fragment from.
     *  \param frag the vertices of G that are part of the fragment.
     *  \param overlap the vertices of G that are part of the overlap region. */
    Fragment(graph::MolecularGraph& G,
             std::vector<graph::MGVertex>& frag,
             std::vector<graph::MGVertex>& overlap);
    /// \todo Add fragment generation directly from molecule.
    ///       Needs shared_ptr impl of mol/atom etc data.
    
    Fragment(const Fragment& frag);
    Fragment(Fragment&& frag);
    Fragment& operator=(const Fragment& frag);
    Fragment& operator=(Fragment&& frag);
    
    graph::CondensedMolecularGraph& GetGraph() const;
    const std::vector<graph::CMGVertex>& GetFragment() const;
    size_t Size() const;
    const std::vector<OverlapVertex>& GetOverlap() const;
    bool IsFragmentVertex(graph::CMGVertex& v) const;
    bool IsOverlapVertex(graph::CMGVertex& v) const;
    const std::vector<AtmType>& GetAtoms() const;
    const std::vector<BndType>& GetBonds() const;
    const std::vector<AngType>& GetAngles() const;
    const std::vector<DhdType>& GetDihedrals() const;
    
    bool operator==(const Fragment& frag) const;
    bool operator!=(const Fragment& frag) const { return !((*this) == frag); }
    bool operator<(const Fragment& frag) const;
    bool operator>(const Fragment& frag) const;
    bool operator<=(const Fragment& frag) const { return !((*this) > frag); }
    bool operator>=(const Fragment& frag) const { return !((*this) < frag); }
    
  private:
    struct FragmentData;
    std::shared_ptr<FragmentData> _dat;
  };
  
  class Athenaeum {
    friend class cereal::access;
  public:
    struct Settings {
      // Maximum vertex count for automatic fragment generation
      static uint32_t AtomLimit;
      static uint32_t MinimumFragmentSize;
      static uint32_t MaximumFragmentSize;
      // Fully fragment cycles
      static bool FragmentCycles;
      // Maximum cycle size to be regarded as a cycle
      static uint32_t MaximumCycleSize;
      static uint32_t DefaultOverlap;
      static uint32_t DefaultCycleOverlap;
    };

    using FragContain = std::vector<Fragment>;
    using MoleculeFragments = std::map<sMolecule, FragContain>;

  private:
    template <class Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    Athenaeum() = default;
    Athenaeum(const Athenaeum& ath);
    Athenaeum(Athenaeum&& ath);
    Athenaeum& operator=(const Athenaeum& ath);
    Athenaeum& operator=(Athenaeum&& ath);
    Athenaeum(Forcefield& ff);
    Athenaeum(Forcefield& ff, uint32_t overlap);
    Athenaeum(Forcefield& ff, uint32_t overlap, uint32_t ring_overlap);

    size_t NumFragments() const;
    size_t NumFragments(Molecule& mol) const;

    const MoleculeFragments& GetFragments() const;
    const FragContain& GetFragments(Molecule& mol) const;
    bool HasFragments(Molecule& mol) const;

    const Forcefield& GetForcefield() const { return _ff; }
    bool IsSelfConsistent() const { return _man; }
    void SetSelfConsistent() { _man = true; }
//    bool CheckSelfConsistent();
    
    /*! \brief Adds the given fragment.
     *  \returns if the fragment was added or not. */
    bool AddFragment(Molecule& mol, Fragment& frag);
    
    /*! \brief Determines all the fragments of a molecule and adds them.
     *  \returns the number of fragments added. */
    size_t AddAllFragments(Molecule& mol);

  private:
    //! \brief Forcefield used
    Forcefield _ff;
    //! \brief Overlap length
    uint32_t _overlap;
    //! \brief Ring overlap length
    uint32_t _roverlap;
    //! \brief Manual
    bool _man;
    //! \brief fragments
    MoleculeFragments _frags;
  };
  
  void SaveAthenaeum(const Athenaeum& ath, std::string path);
  Athenaeum LoadAthenaeum(std::string path);
  
}

#endif /* INDIGOX_CLASSES_ATHENAEUM_HPP */
