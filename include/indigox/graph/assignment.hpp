/*! \file assignment.hpp */
#ifndef INDIGOX_GRAPH_ASSIGNMENT_HPP
#define INDIGOX_GRAPH_ASSIGNMENT_HPP

#include <limits>
#include <map>
#include <memory>
#include <vector>

#include "base_graph.hpp"
#include "../utils/numerics.hpp"

// Foward declares
namespace indigox {
  namespace test { class IXAssignmentGraph; }
  namespace graph {
    class IXMolecularGraph;
    class IXAssignmentGraph;
    class IXMGVertex;
    class IXAGVertex;
    class IXMGEdge;
    
    using MolecularGraph = std::shared_ptr<IXMolecularGraph>;
    //! \brief shared_ptr for normal use of the IXAssignmentGraph class.
    using AssignmentGraph = std::shared_ptr<IXAssignmentGraph>;
    using MGVertex = std::shared_ptr<IXMGVertex>;
    //! \brief shared_ptr for normal use of the IXEAGVertex class.
    using AGVertex = std::shared_ptr<IXAGVertex>;
    using MGEdge = std::shared_ptr<IXMGEdge>;
    
    using _MolecularGraph = std::weak_ptr<IXMolecularGraph>;
    /*! \brief weak_ptr for non-ownership reference to the IXAssignmentGraph class.
     *  \details Intended for internal use only. */
    using _AssignmentGraph = std::weak_ptr<IXAssignmentGraph>;
    using _MGVertex = std::weak_ptr<IXMGVertex>;
    /*! \brief weak_ptr for non-ownership reference to the IXAGVertex class.
     *  \details Intended for internal use only. */
    using _AGVertex = std::weak_ptr<IXAGVertex>;
    using _MGEdge = std::weak_ptr<IXMGEdge>;
  }
}

// Local declarations
namespace indigox::graph {
  
  /*! \brief Class for the vertices of an IXAssignmentGraph. */
  class IXAGVertex : public std::enable_shared_from_this<IXAGVertex> {
  public:
    IXAGVertex() = delete; // no default constructor
    //! \brief Friendship allows IXAssignmentGraph to construct vertices
    friend class IXAssignmentGraph;
    
  private:
    /*! \brief Construct an IXAGVertex from an MGVertex.
     *  \param v the MGVertex to associate with this vertex. */
    IXAGVertex(const MGVertex v)
    : _source_vert(v), _source_edge(), _pre(0), _count(0) { }
    
    /*! \brief Construct an IXAGVertex from an MGEdge.
     *  \param e the MGEdge to associate with this vertex. */
    IXAGVertex(const MGEdge e)
    : _source_vert(), _source_edge(e), _pre(0), _count(0) { }
    
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
    inline uint_ GetPreAssignedCount() const { return _pre; }
    
    /*! \brief Set the number of pre-assigned eletrons.
     *  \details ElectronAssignerAlgorithms will not modify this value.
     *  \param count the number of electrons to pre-assign. */
    inline void SetPreAssignedCount(const uint_ count) { _pre = count; }
    
    /*! \brief Get the total number of assigned electrons.
     *  \details The total number of electrons is the sum of the number of
     *  pre-assigned electrons and the number of assigned electrons.
     *  \return the total number of assigned electrons. */
    inline uint_ GetTotalAssigned() const { return _pre + _count; }
    
    /*! \brief Set the number of assigned eletrons.
     *  \param count the number of electrons to assign. */
    inline void SetAssignedCount(const uint_ count) { _count = count; }
    
  private:
    //! \brief reference to an MGVertex
    _MGVertex _source_vert;
    //! \brief reference to an MDEdge
    _MGEdge _source_edge;
    /*! \property _pre
     *  \brief pre-assigned electrons
     *  \property _count
     *  \brief assigned electrons */
    uint_ _pre, _count;
  };
  
  /*! \brief Class used to assign electrons to a molecule.
   *  \details An IXAssignmentGraph is generated from an IXMolecularGraph at a
   *  fixed point in time. This means that there are no modification methods
   *  available for adding or removing vertices to/from an IXAssignmentGraph.
   *  Additionally, the edges of the graph are not accessible in any way. */
  class IXAssignmentGraph {
    //! \brief Friendship allows IXAssignmentGraph to be properly tested
    friend class indigox::test::IXAssignmentGraph;
    
    //! \brief Type of the internally utilised graph.
    using graph_type = IXGraphBase<IXAGVertex, std::nullptr_t>;
    
  public:
    //! \brief Type of the iterator returned by the GetVertices() method.
    using VertIter = std::vector<AGVertex>::const_iterator;
    //! \brief Type of the iterator returned by the GetNeihbours() method.
    using NbrsIter = std::vector<AGVertex>::const_iterator;
    
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
      return _verts_v.find(v) != _verts_v.end();
    }
    
    /*! \brief Check if an IXMGEdge has an associated vertex in this graph.
     *  \param e the edge to check for.
     *  \return if \p e is associated with the graph or not. */
    inline bool HasVertex(const MGEdge& e) const {
      return _verts_e.find(e) != _verts_e.end();
    }
    
    /*! \brief Get the AGVertex associated with an MGVertex.
     *  \details If the provided MGVertex is not associated with the graph, the
     *  returned AGVertex is empty.
     *  \param v the MGVertex to get associated AGVertex for.
     *  \return the AGVertex associated with the MGVertex. */
    inline AGVertex GetVertex(const MGVertex& v) const {
      return HasVertex(v) ? _verts_v.at(v) : AGVertex();
    }
    
    /*! \brief Get the AGVertex associated with an MGEdge.
     *  \details If the provided MGEdge is not associated with the graph, the
     *  returned AGVertex is empty.
     *  \param e the MGEdge to get associated AGVertex for.
     *  \return the AGVertex associated with the MGEdge. */
    inline AGVertex GetVertex(const MGEdge& e) const {
      return HasVertex(e) ? _verts_e.at(e) : AGVertex();
    }
    
    /*! \brief Get the number of vertices in the graph.
     *  \return the number of vertices in the graph. */
    inline size_ NumVertices() const { return _g.NumVertices(); }
    
    /*! \brief Get the degree of a vertex.
     *  \details If the provided vertex is not part of the graph, the returned
     *  value is std::numeric_limits<size_>::max()
     *  \return the degree of the vertex. */
    inline size_ Degree(const AGVertex& v) const {
      return HasVertex(v) ? _g.Degree(v.get()) : std::numeric_limits<size_>::max();
    }
    
    /*! \brief Get the neighbours of a vertex.
     *  \details If the vertex is not part of the graph, the returned iterator
     *  pair are an empty range.
     *  \param v the vertex to get neighbours of.
     *  \return a pair of iterators across the neighbours of the vertex. */
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(const AGVertex& v) const {
      return HasVertex(v)
        ? std::make_pair(_nbrs.at(v).begin(), _nbrs.at(v).end())
        : std::make_pair(_nbrs.begin()->second.end(), _nbrs.begin()->second.end());
    }
    
    /*! \brief Get the vertices of the graph
     *  \return a pair of iterators across the vertices of the graph. */
    inline std::pair<VertIter, VertIter> GetVertices() const {
      return {_verts.begin(), _verts.end()};
    }
    
    /*! \brief Check if the graph is connected.
     *  \return if the graph is connected or not. */
    inline bool IsConnected() { return _g.NumConnectedComponents() == 1; }
    
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
    std::map<MGVertex, AGVertex> _verts_v;
    //! \brief Map MGEdges to their corresponding AGVertex
    std::map<MGEdge, AGVertex> _verts_e;
    //! \brief All vertices in the graphs
    std::vector<AGVertex> _verts;
    //! \brief Neighbours of vertices
    std::map<AGVertex, std::vector<AGVertex>> _nbrs;
  };
  
}

#endif /* INDIGOX_GRAPH_ASSIGNMENT_HPP */
