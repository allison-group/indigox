//
//  permutablegraph.hpp
//  indigox
//
//  Created by Ivan Welsh on 13/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//

#ifndef INDIGOX_CLASSES_PERMUTABLEGRAPH_HPP
#define INDIGOX_CLASSES_PERMUTABLEGRAPH_HPP

#include <memory>
#include <vector>

#include "molecular_graph.hpp"
#include "../utils/graph.hpp"

namespace indigox {
  struct PermVertProp {
    MolVertPair source;
    uid_t bag;
  };
  
  class _PermutableGraph;
  typedef std::shared_ptr<_PermutableGraph> PermutableGraph;
  typedef utils::Graph<PermVertProp> _PermGraph;
  typedef _PermGraph::VertType PermVertex;
  typedef _PermGraph::NbrsIter PermNbrsIter;
  typedef _PermGraph::VertIter PermVertIter;
  typedef _PermGraph::NbrsIterPair PermNbrsIterPair;
  typedef _PermGraph::VertIterPair PermVertIterPair;
  typedef _PermGraph::EdgeIterPair PermEdgeIterPair;
  
  typedef std::vector<PermVertex> ElimOrder;
  
  class _PermutableGraph : public _PermGraph
  {
  public:
    _PermutableGraph();
    _PermutableGraph(MolecularGraph);
    _PermutableGraph(PermutableGraph);
    
  public:
    void SetInput(MolecularGraph);
    void EliminateVertex(PermVertex);
    std::string ToDGFString();
    MolecularGraph GetSourceGraph() { return source_; }
    std::string PGVToMGVTable();
    
  private:
    MolecularGraph source_;
  };
}

#endif /* INDIGOX_CLASSES_PERMUTABLEGRAPH_HPP */
