//
//  tree_decomposition.hpp
//  indigox
//
//  Created by Welsh, Ivan on 20/11/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_CLASSES_NICETREEDECOMP_HPP
#define INDIGOX_CLASSES_NICETREEDECOMP_HPP

#include "../utils/graph.hpp"
#include "treedecomp.hpp"

#include <set>
#include <string>

namespace indigox {

  struct NTDVertProp {
    uint32_t id;
    std::set<MolVertPair> bag;
    std::pair<char, MolVertPair> kind;
  };

  typedef utils::Graph<NTDVertProp, utils::NoProperty, utils::DirectedGraph>
      _NTDGraph;

  typedef _NTDGraph::VertType NTDVertex;
  typedef _NTDGraph::VertIter NTDVertexIter;
  typedef _NTDGraph::EdgeType NTDEdge;
  typedef _NTDGraph::EdgeIter NTDEdgeIter;
  typedef _NTDGraph::NbrsIter NTDNeighboursIter;
  typedef _NTDGraph::PredIter NTDPredecessorsIter;

  typedef _NTDGraph::VertTypePair NTDVertPair;
  typedef _NTDGraph::VertIterPair NTDVertIterPair;
  typedef _NTDGraph::EdgeIterPair NTDEdgeIterPair;
  typedef _NTDGraph::NbrsIterPair NTDNbrsIterPair;
  typedef _NTDGraph::PredIterPair NTDPredIterPair;
  typedef _NTDGraph::VertBool NTDVertBool;
  typedef _NTDGraph::EdgeBool NTDEdgeBool;

  class _NTDecomp : public _NTDGraph {
  public:
    _NTDecomp();
    _NTDecomp(TDecomp);
    _NTDecomp(TDecomp, TDVertex);

    void SetInput(TDecomp);
    void SetInput(TDecomp, TDVertex);
    void TopologicalSort(std::vector<NTDVertex> &) const;
    TDecomp GetSourceGraph() {
      return source_;
    }
    size_t GetWidth();

  private:
    uint32_t vert_count_ = 0;
    TDecomp source_;
  };

  typedef std::shared_ptr<_NTDecomp> NTDecomp;

} // namespace indigox

#endif /* INDIGOX_CLASSES_NICETREEDECOMP_HPP */
