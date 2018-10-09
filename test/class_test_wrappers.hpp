#ifndef INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP
#define INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/quad.hpp>
#include <indigox/utils/triple.hpp>

namespace indigox::test {
  struct IXAssignmentGraph {
    graph::IXAssignmentGraph g;
  public:
    using VertIter = graph::IXAssignmentGraph::VertIter;
    using NbrsIter = graph::IXAssignmentGraph::NbrsIter;
    // private wrapping members
    IXAssignmentGraph() : g(graph::MolecularGraph()) { }
    inline void AddEdges(graph::MGVertex s, graph::MGVertex t, graph::MGEdge e) { g.AddEdges(s,t,e); }
    inline graph::AGVertex AddVertex(graph::MGVertex v) { return g.AddVertex(v); }
    inline void DetermineAllNeighbours() { g.DetermineAllNeighbours(); }
    
    // public wrapping members
    IXAssignmentGraph(graph::MolecularGraph G) : g(G) { }
    inline size_ Degree(graph::AGVertex v) { return g.Degree(v); }
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(graph::AGVertex v) { return g.GetNeighbours(v); }
    inline graph::AGVertex GetVertex(graph::MGVertex v) { return g.GetVertex(v); }
    inline graph::AGVertex GetVertex(graph::MGEdge v) { return g.GetVertex(v); }
    inline std::pair<VertIter, VertIter> GetVertices() { return g.GetVertices(); }
    inline bool HasVertex(graph::AGVertex v) { return g.HasVertex(v); }
    inline bool HasVertex(graph::MGVertex v) { return g.HasVertex(v); }
    inline bool HasVertex(graph::MGEdge e) { return g.HasVertex(e); }
    inline bool IsConnected() { return g.IsConnected(); }
    inline size_ NumVertices() { return g.NumVertices(); }
  };
  
  struct IXMolecularGraph {
    graph::IXMolecularGraph g;
  public:
    typedef graph::IXMolecularGraph::EdgeIter EdgeIter;
    typedef graph::IXMolecularGraph::NbrsIter NbrsIter;
    typedef graph::IXMolecularGraph::VertIter VertIter;
    typedef graph::IXMolecularGraph::CompIter CompIter;
    
    // private wrapping members
    IXMolecularGraph() : g(Molecule()) { }
    IXMolecularGraph(Molecule m) : g(m) {}
    inline graph::MGEdge AddEdge(Bond bnd) { return g.AddEdge(bnd); }
    inline graph::MGVertex AddVertex(Atom atm) { return g.AddVertex(atm); }
    inline void Clear() { g.Clear(); }
    inline void RemoveEdge(graph::MGEdge e) { g.RemoveEdge(e); }
    inline void RemoveEdge(graph::MGVertex u, graph::MGVertex v) { g.RemoveEdge(u,v); }
    inline void RemoveVertex(graph::MGVertex v) { g.RemoveVertex(v); }
    
    // public wrapping methods
    inline size_ Degree(graph::MGVertex v) { return g.Degree(v); }
    inline std::pair<CompIter, CompIter> GetConnectedComponents() { return g.GetConnectedComponents(); }
    inline graph::MGEdge GetEdge(graph::MGVertex u, graph::MGVertex v) { return g.GetEdge(u, v); }
    inline graph::MGEdge GetEdge(Bond b) { return g.GetEdge(b); }
    inline std::pair<EdgeIter, EdgeIter> GetEdges() { return g.GetEdges(); }
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(graph::MGVertex v) { return g.GetNeighbours(v); }
    inline graph::MGVertex GetSource(graph::MGEdge e) { return g.GetSource(e); }
    inline graph::MGVertex GetTarget(graph::MGEdge e) { return g.GetTarget(e); }
    inline graph::MGVertex GetVertex(Atom a) { return g.GetVertex(a); }
    inline std::pair<graph::MGVertex, graph::MGVertex> GetVertices(graph::MGEdge e) { return g.GetVertices(e); }
    inline std::pair<VertIter, VertIter> GetVertices() { return g.GetVertices(); }
    inline bool HasEdge(Bond b) { return g.HasEdge(b); }
    inline bool HasEdge(graph::MGEdge e) { return g.HasEdge(e); }
    inline bool HasEdge(graph::MGVertex u, graph::MGVertex v) { return g.HasEdge(u,v); }
    inline bool HasVertex(Atom v) { return g.HasVertex(v); }
    inline bool HasVertex(graph::MGVertex v) { return g.HasVertex(v); }
    inline bool IsConnected() { return g.IsConnected(); }
    inline size_ NumConnectedComponents() { return g.NumConnectedComponents(); }
    inline size_ NumEdges() { return g.NumEdges(); }
    inline size_ NumVertices() { return g.NumVertices(); }
  };
  
  
}

#endif /* INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP */
