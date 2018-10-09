/*! \file assignment.hpp */
#ifndef INDIGOX_GRAPH_ASSIGNMENT_HPP
#define INDIGOX_GRAPH_ASSIGNMENT_HPP
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <vector>

#include <EASTL/vector_map.h>

#include "base_graph.hpp"
#include "../utils/common.hpp"
#include "../utils/fwd_declares.hpp"

// Local declarations
namespace indigox::graph {
  
  /*! \brief Class for the vertices of an IXAssignmentGraph. */
  class IXAGVertex : public std::enable_shared_from_this<IXAGVertex> {
  public:
    IXAGVertex() = delete; // no default constructor
    //! \brief Friendship allows IXAssignmentGraph to construct vertices
    friend class IXAssignmentGraph;
    //! \brief Friendship allows serialisation
    
  private:
    /*! \brief Construct an IXAGVertex from an MGVertex.
     *  \param v the MGVertex to associate with this vertex. */
    IXAGVertex(const MGVertex v)
    : _source_vert(v), _source_edge(), _pre(0), _count(0) { }
    
    /*! \brief Construct an IXAGVertex from an MGEdge.
     *  \details As all bonds must have an order of at least one, the
     *  preassigned count of all edge mapped vertices is always set to two.
     *  \param e the MGEdge to associate with this vertex. */
    IXAGVertex(const MGEdge e)
    : _source_vert(), _source_edge(e), _pre(2), _count(0) { }
    
  public:
    /*! \brief Is the associated IXMolecularGraph member a vertex.
     *  \return if the \p _source_vert member is not empty. */
    inline bool IsVertexMapped() const { return !_source_vert.expired(); }
    
    /*! \brief Is the associated IXMoleculerGraph member an edge.
     *  \return if the \p _source_edge member is not empty. */
    inline bool IsEdgeMapped() const { return !_source_edge.expired(); }
    
    /*! \brief Obtain the referenced IXMolecularGraph vertex.
     *  \details If the IXAGVertex instance is not mapped to a vertex, the
     *  returned MGVertex will be empty.
     *  \return the MGVertex instance associated with this IXAGVertex. */
    inline MGVertex GetSourceVertex() const {
      return IsVertexMapped() ? _source_vert.lock() : MGVertex();
    }
    
    /*! \brief Obtain the referenced IXMolecularGraph edge.
     *  \details If the IXAGVertex instance is not mapped to an edge, the
     *  returned MGVertex will be empty.
     *  \return the MGEdge instance associated with this IXAGVertex. */
    inline MGEdge GetSourceEdge() const {
      return IsEdgeMapped() ? _source_edge.lock() : MGEdge();
    }
    
    /*! \brief Get the number of pre-assigned electrons.
     *  \return the number of pre-assigned eletrons. */
    inline uint32_t GetPreAssignedCount() const { return _pre; }
    
    /*! \brief Set the number of pre-assigned eletrons.
     *  \details ElectronAssignerAlgorithms will not modify this value.
     *  \param count the number of electrons to pre-assign. */
    inline void SetPreAssignedCount(const uint32_t count) { _pre = count; }
    
    /*! \brief Get the number of assigned electrons.
     *  \details Does not include preassigned electrons.
     *  \return the number of assigned electrons. */
    inline uint32_t GetAssignedCount() const { return _count; }
    
    /*! \brief Get the total number of assigned electrons.
     *  \details The total number of electrons is the sum of the number of
     *  pre-assigned electrons and the number of assigned electrons.
     *  \return the total number of assigned electrons. */
    inline uint32_t GetTotalAssigned() const { return _pre + _count; }
    
    /*! \brief Set the number of assigned eletrons.
     *  \param count the number of electrons to assign. */
    inline void SetAssignedCount(const uint32_t count) { _count = count; }
    
  private:
    //! \brief reference to an MGVertex
    _MGVertex _source_vert;
    //! \brief reference to an MDEdge
    _MGEdge _source_edge;
    /*! \property _pre
     *  \brief pre-assigned electrons
     *  \property _count
     *  \brief assigned electrons */
    uint32_t _pre, _count;
  };
  
  /*! \brief Class used to assign electrons to a molecule.
   *  \details An IXAssignmentGraph is generated from an IXMolecularGraph at a
   *  fixed point in time. This means that there are no modification methods
   *  available for adding or removing vertices to/from an IXAssignmentGraph.
   *  Additionally, the edges of the graph are not accessible in any way. */
  class IXAssignmentGraph {
    //! \brief Friendship allows IXAssignmentGraph to be properly tested
    friend struct indigox::test::TestAssignmentGraph;
    
    //! \brief Type of the internally utilised graph.
    using graph_type = IXGraphBase<IXAGVertex, std::nullptr_t>;
    //! \brief Container for vertices
    using VertContain = std::vector<AGVertex>;
    
  public:
    //! \brief Type of the iterator returned by the GetVertices() method.
    using VertIter = VertContain::const_iterator;
    //! \brief Type of the iterator returned by the GetNeihbours() method.
    using NbrsIter = VertContain::const_iterator;
    
  public:
    IXAssignmentGraph() = delete; // no default constructor
    
  public:
    /*! \brief Construct from a molecular graph.
     *  \param g molecular graph to construct from. */
    IXAssignmentGraph(MolecularGraph g);
    
    /*! \brief Check if an IXAGVertex belongs to this graph.
     *  \param v the vertex to check for.
     *  \return if \p v is part of the graph or not. */
    inline bool HasVertex(const AGVertex& v) const {
      return _g.HasVertex(v.get());
    }
    
    /*! \brief Check if an IXMGVertex has an associated vertex in this graph.
     *  \param v the IXMGVertex to check for.
     *  \return if \p v is associated with the graph or not. */
    inline bool HasVertex(const MGVertex& v) const {
      return _v2v.find(v) != _v2v.end();
    }
    
    /*! \brief Check if an IXMGEdge has an associated vertex in this graph.
     *  \param e the edge to check for.
     *  \return if \p e is associated with the graph or not. */
    inline bool HasVertex(const MGEdge& e) const {
      return _e2v.find(e) != _e2v.end();
    }
    
    /*! \brief Get the AGVertex associated with an MGVertex.
     *  \details If the provided MGVertex is not associated with the graph, the
     *  returned AGVertex is empty.
     *  \param v the MGVertex to get associated AGVertex for.
     *  \return the AGVertex associated with the MGVertex. */
    inline AGVertex GetVertex(const MGVertex& v) const {
      return HasVertex(v) ? _v2v.at(v) : AGVertex();
    }
    
    /*! \brief Get the AGVertex associated with an MGEdge.
     *  \details If the provided MGEdge is not associated with the graph, the
     *  returned AGVertex is empty.
     *  \param e the MGEdge to get associated AGVertex for.
     *  \return the AGVertex associated with the MGEdge. */
    inline AGVertex GetVertex(const MGEdge& e) const {
      return HasVertex(e) ? _e2v.at(e) : AGVertex();
    }
    
    /*! \brief Get the number of vertices in the graph.
     *  \return the number of vertices in the graph. */
    inline size_t NumVertices() const { return _g.NumVertices(); }
    
    /*! \brief Get the degree of a vertex.
     *  \details If the provided vertex is not part of the graph, the returned
     *  value is std::numeric_limits<size_>::max()
     *  \return the degree of the vertex. */
    inline size_t Degree(const AGVertex& v) const {
      return HasVertex(v) ? _g.Degree(v.get()) : std::numeric_limits<size_t>::max();
    }
    
    /*! \brief Get the neighbours of a vertex.
     *  \details If the vertex is not part of the graph, the returned iterator
     *  pair are an empty range.
     *  \param v the vertex to get neighbours of.
     *  \return a pair of iterators across the neighbours of the vertex. */
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(const AGVertex& v) const {
      return HasVertex(v) ? std::make_pair(_n.at(v).begin(), _n.at(v).end())
        : std::make_pair(_v.end(), _v.end());
    }
    
    /*! \brief Get the vertices of the graph
     *  \return a pair of iterators across the vertices of the graph. */
    inline std::pair<VertIter, VertIter> GetVertices() const {
      return {_v.begin(), _v.end()};
    }
    
    /*! \brief Check if the graph is connected.
     *  \return if the graph is connected or not. */
//    inline bool IsConnected() { return _g.NumConnectedComponents() == 1; }
    
    /*! \brief Populate the pre-assigned state of all vertices.
     *  \details The assigned electrons are as follows:
     *
     *  - Six electrons on one coordinate F, Cl, and Br.
     *  - Four electrons on one coordinate O, and S.
     *  - Two electrons on one coordinate N.
     *  - Four electrons on two coordinate O, and S.
     *  - No electrons on all other types of atoms.
     *  - Two electrons on all bonds. */
    void PreassignElectrons();
    
  private:
    /*! \brief Add a vertex to the graph.
     *  \param v the MGVertex associated with the new vertex.
     *  \return the newly added vertex. */
    AGVertex AddVertex(MGVertex v);
    
    /*! \brief Add edges to the graph.
     *  \details Two edges are added. The first from \p s to \p e and the second
     *  from \p e to \p t. If either \p s or \p t are not already added to the
     *  graph, they will be added.
     *  \param s,t the source and target vertices of \p e.
     *  \param e the edge. */
    void AddEdges(MGVertex s, MGVertex t, MGEdge e);
    
    //! \brief Populate the _nbrs member with all the correct neighbours.
    void DetermineAllNeighbours();
    
  private:
    //! \brief Source molecular graph
    _MolecularGraph _source;
    //! \brief Underlying graph
    graph_type _g;
    //! \brief Map MGVertices to their corresponding AGVertex
    eastl::vector_map<MGVertex, AGVertex> _v2v;
    //! \brief Map MGEdges to their corresponding AGVertex
    eastl::vector_map<MGEdge, AGVertex> _e2v;
    //! \brief All vertices in the graphs
    VertContain _v;
    //! \brief Neighbours of vertices
    std::map<AGVertex, VertContain> _n;
  };
  
}

#endif /* INDIGOX_GRAPH_ASSIGNMENT_HPP */
