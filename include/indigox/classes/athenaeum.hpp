#include <map>
#include <vector>


#include "molecule.hpp"
#include "../graph/condensed.hpp"

#ifndef INDIGOX_CLASSES_ATHENAEUM_HPP
#define INDIGOX_CLASSES_ATHENAEUM_HPP

namespace indigox {
  
  // Forward declares
  class IXFragment;
  using Fragment = std::shared_ptr<IXFragment>;
  
  class IXAthenaeum;
  using Athenaeum = std::shared_ptr<IXAthenaeum>;
  
  class IXFragment {
  public:
//    /*! \brief Type of overlapping vertex.
//     *  \details Each vertex of an overlap is given a type. This is to enable
//     *  more detailed isomorphism matching methods to be used when matching. */
//    enum class OverlapType {
//      Distance_1,
//      Distance_2,
//      Distance_3,
//      Distance_4,
//      Distance_5,
//      Distance_6,
//      Distance_7,
//      Distance_8,
//      Distance_9,
//      RingDistance_1,
//      RingDistance_2,
//      RingDistance_3,
//      RingDistance_4,
//      RingDistance_5,
//      RingDistance_6,
//      RingDistance_7,
//      RingDistance_8,
//      RingDistance_9
//    };
    
    //! \brief Type used to match an overlap vertex to its type of overlap.
    using OverlapVertex = graph::CMGVertex; //std::pair<OverlapType, graph::CMGVertex>;
    friend bool operator==(const IXFragment& a, const IXFragment& b);
    
  public:
    IXFragment() = delete;  // no default constructor
    
    /*! \brief Normal constructor
     *  \details Constructs a fragment of the given graph from the combined
     *  fragment and overlap vertices.
     *  \param G the graph to generate fragment from.
     *  \param frag the vertices of G that are part of the fragment.
     *  \param overlap the vertices of G that are part of the overlap region. */
    template <class Input>
    IXFragment(const graph::CondensedMolecularGraph& G,
               const Input& frag, const Input& overlap)
    : _frag(frag.begin(), frag.end()), _overlap(overlap.begin(), overlap.end()) {
      std::vector<graph::CMGVertex> tmp(frag.begin(), frag.end());
      tmp.insert(tmp.end(), overlap.begin(), overlap.end());
      _g = G->InduceSubgraph(tmp.begin(), tmp.end());
    }
    
  private:
    //! \brief Fragment graph, including overlap
    graph::CondensedMolecularGraph _g;
    //! \brief Vertices of the fragment
    std::vector<graph::CMGVertex> _frag;
    //! \brief Vertices of the overlap
    std::vector<OverlapVertex> _overlap;
  };
  
  inline bool operator==(const IXFragment& a, const IXFragment& b) {
    return (a._frag.size() == b._frag.size()
            && a._overlap.size() == b._overlap.size()
            && a._frag == b._frag
            && a._overlap == b._overlap);
  }
  inline bool operator==(const Fragment& a, const Fragment& b) { return *a == *b; }
  
  class IXAthenaeum {
  public:
    struct Settings {
      // Maximum vertex count for automatic fragment generation
      static uint_ AutomaticVertexLimit;
      // Fully fragment cycles
      static bool AutomaticFragmentCycles;
      // Maximum cycle size to be regarded as a cycle
      static uint_ AutomaticMaximumCycleSize;
      
    };
    
    using MolFrags = std::vector<Fragment>;
    using CondensedGraphs = std::map<Molecule, graph::CondensedMolecularGraph>;
    using FragStore = std::map<graph::CondensedMolecularGraph, MolFrags>;
    
    IXAthenaeum();
    IXAthenaeum(uint_ overlap, uint_ ring_overlap);
    
    /*! \brief Determines all the fragments of a molecule and adds them.
     *  \returns the number of fragments added. */
    size_ AddAllFragments(const Molecule& mol);
    /*! \brief Adds the given fragment.
     *  \returns if the fragment was added or not. */
    bool AddFragment(const Molecule& mol, Fragment frag);
    
    size_ NumFragments() const;
    size_ NumFragments(const Molecule& mol) const;
    
    const MolFrags& GetFragments(const Molecule& mol) const;
    const FragStore& GetFragments() const { return _frags; }
    
    bool HasFragments(const Molecule& mol) const {
      return _graphs.find(mol) != _graphs.end();
    }
    
  private:
    //! \brief Overlap length
    uint_ _overlap;
    //! \brief Ring overlap length
    uint_ _roverlap;
    //! \brief Condensed graphs
    CondensedGraphs _graphs;
    //! \brief Per molecule fragments
    FragStore _frags;
  };
  
}

#endif /* INDIGOX_CLASSES_ATHENAEUM_HPP */
