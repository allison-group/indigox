#include <algorithm>
#include <map>
#include <vector>

#include "../utils/fwd_declares.hpp"
#include "../utils/triple.hpp"
#include "../utils/quad.hpp"
#include "molecule.hpp"
#include "../graph/condensed.hpp"

#ifndef INDIGOX_CLASSES_ATHENAEUM_HPP
#define INDIGOX_CLASSES_ATHENAEUM_HPP

namespace indigox {
  
  class Fragment {
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
    
    friend bool operator==(const Fragment& a, const Fragment& b);
    
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
    
    Fragment(const Fragment& frag);
    Fragment(Fragment&& frag);
    Fragment& operator=(const Fragment& frag);
    Fragment& operator=(Fragment&& frag);
    
    graph::CondensedMolecularGraph& GetGraph() const;
    const std::vector<graph::CMGVertex>& GetFragment() const;
    const std::vector<OverlapVertex>& GetOverlap() const;
    bool IsFragmentVertex(graph::CMGVertex& v) const;
    bool IsOverlapVertex(graph::CMGVertex& v) const;
    const std::vector<AtmType>& GetAtoms() const;
    const std::vector<BndType>& GetBonds() const;
    const std::vector<AngType>& GetAngles() const;
    const std::vector<DhdType>& GetDihedrals() const;
    
  private:
    struct FragmentData;
    std::shared_ptr<FragmentData> _dat;
  };
  
//  class Athenaeum {
//  public:
//    struct Settings {
//      // Maximum vertex count for automatic fragment generation
//      static uint32_t AutomaticVertexLimit;
//      // Fully fragment cycles
//      static bool AutomaticFragmentCycles;
//      // Maximum cycle size to be regarded as a cycle
//      static uint32_t AutomaticMaximumCycleSize;
//      
//    };
//    
//    using MolFrags = std::vector<Fragment>;
//    using CondensedGraphs = std::map<Molecule, graph::CondensedMolecularGraph>;
//    using FragStore = std::map<graph::CondensedMolecularGraph, MolFrags>;
//    
//    Athenaeum() = delete;
//    Athenaeum(Forcefield& ff);
//    Athenaeum(Forcefield& ff, uint32_t overlap, uint32_t ring_overlap);
//    
//    /*! \brief Determines all the fragments of a molecule and adds them.
//     *  \returns the number of fragments added. */
//    size_t AddAllFragments(Molecule& mol);
//    
//    /*! \brief Adds the given fragment.
//     *  \returns if the fragment was added or not. */
//    bool AddFragment(Molecule& mol, Fragment frag);
//    
//    size_t NumFragments() const;
//    size_t NumFragments(Molecule& mol) const;
//    
//    const MolFrags& GetFragments(Molecule& mol) const;
//    const FragStore& GetFragments() const { return _frags; }
//    
//    bool HasFragments(Molecule& mol) const {
//      return _graphs.find(mol) != _graphs.end();
//    }
//    
//    Forcefield GetForcefield() const;
//    
//    bool IsManualSelfConsistent() const { return _man; }
//    
//    void SetManual() { _man = true; }
//    
//  private:
//    //! \brief Forcefield used
//    sForcefield _ff;
//    //! \brief Overlap length
//    uint32_t _overlap;
//    //! \brief Ring overlap length
//    uint32_t _roverlap;
//    //! \brief Manual
//    bool _man;
//    //! \brief Condensed graphs
//    CondensedGraphs _graphs;
//    //! \brief Per molecule fragments
//    FragStore _frags;
//  };
  
}

#endif /* INDIGOX_CLASSES_ATHENAEUM_HPP */
