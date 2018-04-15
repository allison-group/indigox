//
//  molecular_graph.hpp
//  indigox
//
//  Created by Welsh, Ivan on 25/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_CLASSES_MOLECULAR_GRAPH_HPP
#define INDIGOX_CLASSES_MOLECULAR_GRAPH_HPP

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include <boost/logic/tribool.hpp>

#include "../utils/graph.hpp"

namespace indigox {
  
  class IXAtom;
  class IXBond;
  class IXAngle;
  class IXDihedral;
  class IXMolecule;
  typedef std::shared_ptr<IXAtom> Atom;
  typedef std::shared_ptr<IXBond> Bond;
  typedef std::shared_ptr<IXAngle> Angle;
  typedef std::shared_ptr<IXDihedral> Dihedral;
  typedef std::shared_ptr<IXMolecule> Molecule;
  typedef std::weak_ptr<IXAtom> _Atom;
  typedef std::weak_ptr<IXBond> _Bond;
  typedef std::weak_ptr<IXAngle> _Angle;
  typedef std::weak_ptr<IXDihedral> _Dihedral;
  typedef std::weak_ptr<IXMolecule> _Molecule;
  
  struct MolVertProp {
    Atom atom;
    int component = -1;
  };
  
  struct MolEdgeProp {
    Bond bond;
  };
  
  // TODO: Replace OutEdgeList and VertexList with boost::vecS. Doing so
  //       will require handling of all invalidated descriptors on removal
  typedef utils::Graph<MolVertProp, MolEdgeProp> _MolGraph;
  typedef _MolGraph::VertType MolVertex;
  typedef _MolGraph::VertIter MolVertexIter;
  typedef _MolGraph::EdgeType MolEdge;
  typedef _MolGraph::EdgeIter MolEdgeIter;
  typedef _MolGraph::NbrsIter MolNeighboursIter;
  
  typedef _MolGraph::VertTypePair MolVertPair;
  typedef _MolGraph::VertIterPair MolVertIterPair;
  typedef _MolGraph::EdgeIterPair MolEdgeIterPair;
  typedef _MolGraph::NbrsIterPair MolNbrsIterPair;
  typedef _MolGraph::VertBool MolVertBool;
  typedef _MolGraph::EdgeBool MolEdgeBool;
  typedef std::tuple<MolVertex, MolVertex, MolEdgeProp> MolVertPairEdgeProp;
  
  typedef std::map<MolVertex, int> MolVertIdxMap;
  typedef std::map<MolEdge, int> MolEdgeIdxMap;
  
  class _MolecularGraph : public _MolGraph {
  private:
    boost::tribool planar_ = boost::indeterminate;
    int16_t totalCharge_ = 0;
    size_t num_components_ = 1;
    Molecule source_;
    
  public:
    using _MolGraph::AddVertex;
    using _MolGraph::GetEdge;
    using _MolGraph::AddEdge;
    
    _MolecularGraph();
    _MolecularGraph(Molecule);
    
    MolVertex AddVertex(Atom);
    MolEdgeBool AddEdge(MolVertex, MolVertex, Bond);
    
    // Graph operations
    inline int16_t GetTotalCharge() const { return totalCharge_; }
    inline void SetTotalCharge(int16_t q) { totalCharge_ = q; }
    
    /// @todo transfer connected components to underlying graph
    size_t NumConnectedComponents();
    
    std::string ToDGFString();
    
  };
  typedef std::shared_ptr<_MolecularGraph> MolecularGraph;
  
}  // namespace indigox

#endif /* INDIGOX_CLASSES_MOLECULAR_GRAPH_HPP */
