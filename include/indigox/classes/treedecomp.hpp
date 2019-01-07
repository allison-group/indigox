//
//  tree_decomposition.hpp
//  indigox
//
//  Created by Welsh, Ivan on 11/01/18.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_CLASSES_TREEDECOMP_HPP
#define INDIGOX_CLASSES_TREEDECOMP_HPP

#include "../utils/graph.hpp"
#include "permutablegraph.hpp"

#include <map>
#include <set>
#include <string>

namespace indigox {

  struct TDVertProp {
    std::set<PermVertex> bag;
  };

  typedef utils::Graph<TDVertProp, utils::NoProperty, utils::UndirectedGraph>
      _TDGraph;

  typedef _TDGraph::VertType TDVertex;
  typedef _TDGraph::VertIter TDVertexIter;
  typedef _TDGraph::EdgeType TDEdge;
  typedef _TDGraph::EdgeIter TDEdgeIter;
  typedef _TDGraph::NbrsIter TDNeighboursIter;
  typedef _TDGraph::PredIter TDPredecessorsIter;

  typedef _TDGraph::VertTypePair TDVertPair;
  typedef _TDGraph::VertIterPair TDVertIterPair;
  typedef _TDGraph::EdgeIterPair TDEdgeIterPair;
  typedef _TDGraph::NbrsIterPair TDNbrsIterPair;
  typedef _TDGraph::PredIterPair TDPredIterPair;
  typedef _TDGraph::VertBool TDVertBool;
  typedef _TDGraph::EdgeBool TDEdgeBool;

  //  class _TDecomp

  class _TDecomp : public _TDGraph {
  public:
    _TDecomp();
    _TDecomp(PermutableGraph, ElimOrder &);

    void SetInput(PermutableGraph, ElimOrder &);
    std::string ToString();
    size_t GetWidth() const {
      return upperBound_;
    }
    inline PermutableGraph GetSourceGraph() {
      return originalG_;
    }

  private:
    void FromOrder(ElimOrder &, size_t);

  private:
    PermutableGraph permG_, originalG_;
    std::map<PermVertex, TDVertex> bagMap_;
    std::map<PermVertex, uid_t> elimiIndex_;
    size_t upperBound_;
  };
  typedef std::shared_ptr<_TDecomp> TDecomp;

} // namespace indigox

#endif /* INDIGOX_CLASSES_TREEDECOMP_HPP */
