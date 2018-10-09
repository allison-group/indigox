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
    using Vert = graph::MGVertex;
    friend bool operator==(const IXFragment& a, const IXFragment& b);
    
  public:
    IXFragment() = delete;  // no default constructor
    
    /*! \brief Normal constructor
     *  \details Constructs a fragment of the given graph from the combined
     *  fragment and overlap vertices.
     *  \param G the graph to generate fragment from.
     *  \param frag the vertices of G that are part of the fragment.
     *  \param overlap the vertices of G that are part of the overlap region. */
    IXFragment(const graph::CondensedMolecularGraph& G,
               const Molecule& mol,
               const std::vector<graph::CMGVertex>& frag,
               const std::vector<graph::CMGVertex>& overlap);
    
    graph::CondensedMolecularGraph GetGraph() const { return _g; }
    
    const std::vector<graph::CMGVertex>& GetFragmentVertices() const {
      return _frag; }
    
    const std::vector<OverlapVertex>& GetOverlapVertices() const {
      return _overlap; }
    
    bool IsFragmentVertex(graph::CMGVertex& v) const {
      return std::find(_frag.begin(), _frag.end(), v) != _frag.end();
    }
    
    bool IsOverlapVertex(graph::CMGVertex& v) const {
      return std::find(_overlap.begin(), _overlap.end(), v) != _overlap.end();
    }
    
    const std::vector<Vert>& GetAtomVertices() const { return _atms; }
    const std::vector<std::pair<Vert,Vert>>& GetBondVertices() const {
      return _bnds; }
    const std::vector<stdx::triple<Vert,Vert,Vert>>& GetAngleVertices() const {
      return _angs; }
    const std::vector<stdx::quad<Vert,Vert,Vert,Vert>>&
    GetDihedralVertices() const { return _dhds; }
    
    Molecule GetMolecule() const { return _mol.lock(); }
    
  private:
    //! \brief Fragment graph, including overlap
    graph::CondensedMolecularGraph _g;
    //! \brief Molecule fragment is from
    _Molecule _mol;
    //! \brief Atom vertices
    std::vector<Vert> _atms;
    //! \brief Bond vertices
    std::vector<std::pair<Vert, Vert>> _bnds;
    //! \brief Angle vertices
    std::vector<stdx::triple<Vert, Vert, Vert>> _angs;
    //! \brief Dihedral vertices
    std::vector<stdx::quad<Vert, Vert, Vert, Vert>> _dhds;
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
      static uint32_t AutomaticVertexLimit;
      // Fully fragment cycles
      static bool AutomaticFragmentCycles;
      // Maximum cycle size to be regarded as a cycle
      static uint32_t AutomaticMaximumCycleSize;
      
    };
    
    using MolFrags = std::vector<Fragment>;
    using CondensedGraphs = std::map<Molecule, graph::CondensedMolecularGraph>;
    using FragStore = std::map<graph::CondensedMolecularGraph, MolFrags>;
    
    IXAthenaeum() = delete;
    IXAthenaeum(Forcefield ff);
    IXAthenaeum(Forcefield ff, uint32_t overlap, uint32_t ring_overlap);
    
    /*! \brief Determines all the fragments of a molecule and adds them.
     *  \returns the number of fragments added. */
    size_t AddAllFragments(const Molecule& mol);
    /*! \brief Adds the given fragment.
     *  \returns if the fragment was added or not. */
    bool AddFragment(const Molecule& mol, Fragment frag);
    
    size_t NumFragments() const;
    size_t NumFragments(const Molecule& mol) const;
    
    const MolFrags& GetFragments(const Molecule& mol) const;
    const FragStore& GetFragments() const { return _frags; }
    
    bool HasFragments(const Molecule& mol) const {
      return _graphs.find(mol) != _graphs.end();
    }
    
    Forcefield GetForcefield() const { return _ff; }
    
    bool IsManualSelfConsistent() const { return _man; }
    
    void SetManual() { _man = true; }
    
  private:
    //! \brief Forcefield used
    Forcefield _ff;
    //! \brief Overlap length
    uint32_t _overlap;
    //! \brief Ring overlap length
    uint32_t _roverlap;
    //! \brief Manual
    bool _man;
    //! \brief Condensed graphs
    CondensedGraphs _graphs;
    //! \brief Per molecule fragments
    FragStore _frags;
  };
  
}

#endif /* INDIGOX_CLASSES_ATHENAEUM_HPP */
