/*! \file base.hpp */
#ifndef INDIGOX_GRAPH_BASE_HPP
#define INDIGOX_GRAPH_BASE_HPP

#include <map>
#include <stdexcept>
#include <vector>

#include <boost/bimap.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <indigox/utils/numerics.hpp>

namespace indigox::graph {
  
  //! \brief Type for specifying that a graph is directed.
  struct Directed {
    //! \brief Underlying boost type of a directed graph.
    typedef boost::bidirectionalS is_directed_t;
    //! \brief Boolean that the type is directed.
    static constexpr bool is_directed = true;
  };
  
  //! \brief Type for specifying that a graph is undirected.
  struct Undirected {
    //! \brief Underlying boost type of an undirected graph.
    typedef boost::undirectedS is_directed_t;
    //! \brief Boolean that the type is not directed.
    static constexpr bool is_directed = false;
  };
  
  /*! \brief Template base class for all graphs used in the indigoX library.
   *  \tparam V type of the graph vertices.
   *  \tparam E type of the graph edges.
   *  \tparam D type indicating the directed nature of the graph. Defaults to
   *  undirected. */
  template <class V, class E, class D=Undirected>
  class GraphBase {
  private:
    
    //! \brief Type for applying a label to a vertex or edge.
    struct Label{
      //! \brief An integer label.
      int_ label;
    };
    
    //! \brief Type of the underlying boost graph.
    typedef boost::adjacency_list<boost::setS,               // Edge container
                                  boost::listS,              // Vertex container
                                  typename D::is_directed_t, // Directed nature
                                  Label,                     // Vertex Properties
                                  Label                      // Edge Properties
                                  > G;
    
    // May need to replace G:: with boost::graph_traits<G>::
    //! \brief Type of the graph vertex descriptor.
    typedef typename G::vertex_descriptor VertType;
    //! \brief Type for iterator over graph vertex descriptors.
    typedef typename G::vertex_iterator VertIter;
    //! \brief Type for iterator over neighbours of vertex descriptor.
    typedef typename G::adjacency_iterator NbrsIter;
    //! \brief Type for iterator over predecessors of a vertex descriptor.
    typedef typename G::inv_adjacency_iterator PredIter;
    //! \brief Type of the graph edge descriptor.
    typedef typename G::edge_descriptor EdgeType;
    //! \brief Type for iterator over edges.
    typedef typename G::edge_iterator EdgeIter;
    
    //! \brief Type for bidirectional mapping of V to vertex descriptor type.
    typedef boost::bimap<V, VertType> VertMap;
    //! \brief Type for bidirectional mapping of E to edge descriptor type.
    typedef boost::bimap<E, EdgeType> EdgeMap;
    
  protected:
    //! \brief Underlying boost graph.
    std::shared_ptr<G> _graph;
    //! \brief Map vertices to their descriptors.
    VertMap _verts;
    //! \brief Map edges to their descriptors.
    EdgeMap _edges;
    
  public:
    //! \brief Default constructor
    GraphBase() : _graph(std::make_shared<G>()) { }
    
    /*! \brief Add a new vertex to the graph.
     *  \param v vertex to add.
     *  \throw std::invalid_argument if the vertex already exists within the
     *  graph. */
    void AddVertex(V v) {
      if (HasVertex(v))
        throw std::invalid_argument("Duplicate vertices are not allowed.");
      _verts.insert({v, boost::add_vertex(Label(), *_graph)});
    }
    
    /*! \brief Remove a vertex from the graph.
     *  \details Removing a vertex also removes all edges incident on it.
     *  \param v the vertex to remove.
     *  \throw std::invalid_argument if the vertex does not exist within the
     *  graph. */
    void RemoveVertex(V v) {
      if (!HasVertex(v))
        throw std::invalid_argument("Vertex does not exist.");
      VertType v_ = GetVertexDescriptor(v);
      // Remove adjacent edges
      NbrsIter vi, vi_end;
      std::tie(vi, vi_end) = boost::adjacent_vertices(v_, *_graph);
      for (; vi != vi_end; ++vi) {
        _edges.right.erase(boost::edge(v_, *vi, *_graph).first);
      }
      // Remove incident edges of directed graphs
      if (D::is_directed) {
        PredIter vp, vp_end;
        std::tie(vp, vp_end) = boost::inv_adjacent_vertices(v_, *_graph);
        for (; vp != vp_end; ++vp) {
          _edges.right.erase(boost::edge(*vp, v_, *_graph).first);
        }
      }
      // Remove the vertex
      _verts.right.erase(v_);
      boost::clear_vertex(v_, *_graph);
      boost::remove_vertex(v_, *_graph);
    }
    
    /*! \brief Add a new edge to the graph.
     *  \details Vertex u is used as the source and vertex v as the target. If
     *  u and/or v are not already part of the graph, they are added.
     *  \param u, v vertices the edge is between.
     *  \param e edge to add.
     *  \throw std::invalid_argument if u and v are the same vertex.
     *  \throw std::invalid_argument if there is already an edge between u and
     *  v. */
    void AddEdge(V u, V v, E e) {
      if (u == v)
        throw std::invalid_argument("Adding a self loop is not allowed.");
      if (!HasVertex(u)) AddVertex(u);
      if (!HasVertex(v)) AddVertex(v);
      VertType u_ = GetVertexDescriptor(u);
      VertType v_ = GetVertexDescriptor(v);
      
      std::pair<EdgeType, bool> e_ = boost::edge(u_, v_, *_graph);
      if (e_.second)
        throw std::invalid_argument("Adding parallel edges is not allowed.");
      _edges.insert({e, boost::add_edge(u_, v_, Label(), *_graph).first});
    }
    
    /*! \brief Remove an edge from the graph.
     *  \param e the edge to remove.
     *  \throw std::invalid_argument if the edge does not exist within the
     *  graph. */
    void RemoveEdge(E e) {
      if (!HasEdge(e))
        throw std::invalid_argument("Edge does not exist.");
      EdgeType e_ = GetEdgeDescriptor(e);
      _edges.right.erase(e_);
      boost::remove_edge(e_, *_graph);
    }
    
    /*! \brief Remove an edge from the graph.
     *  \param u, v vertices to remove an edge from between.
     *  \throw std::invalid_argument if no edge between the vertices or at
     *  least one of the vertices does not exist within the graph. */
    void RemoveEdge(V u, V v) {
      if (!HasEdge(u, v))
        throw std::invalid_argument("No such edge exists.");
      VertType u_ = GetVertexDescriptor(u);
      VertType v_ = GetVertexDescriptor(v);
      EdgeType e = boost::edge(u_, v_, *_graph).first;
      _edges.right.erase(e);
      boost::remove_edge(e, *_graph);
    }
    
    /*! \brief Is the vertex in the graph.
     *  \param v vertex to search for.
     *  \return if the requested vertex is contained in the graph or not. */
    bool HasVertex(V v) const {
      return _verts.left.find(v) != _verts.left.end();
    }
    
    /*! \brief Is the edge in the graph.
     *  \param e edge to search for.
     *  \return if the requested edge is contained in the graph or not. */
    bool HasEdge(E e) const { return _edges.left.find(e) != _edges.left.end(); }
    
    /*! \brief Does an edge exist between two vertices.
     *  \param u, v vertices to check between.
     *  \return if there is an edge between the two vertices. */
    bool HasEdge(V u, V v) const {
      if (!HasVertex(u) || !HasVertex(v)) return false;
      VertType u_ = GetVertexDescriptor(u);
      VertType v_ = GetVertexDescriptor(v);
      return boost::edge(u_, v_, *_graph).second;
    }
    
    /*! \brief Number of vertices in the graph.
     *  \return the number of vertices in the graph. */
    size_ NumVertices() const { return boost::num_vertices(*_graph); }
    
    /*! \brief Number of edges in the graph.
     *  \return the number of edges in the graph. */
    size_ NumEdges() const { return boost::num_edges(*_graph); }
    
    //! \brief Removes all edges and vertices from the graph.
    void Clear() {
      _graph->clear();
      _verts.clear();
      _edges.clear();
    }
    
    /*! \brief Degree of a vertex.
     *  \details In the case of a directed graph, the degree of a vertex is the
     *  number of edges leaving the vertex.
     *  \param v the vertex to get the degree of.
     *  \return the degree of the vertex.
     *  \throw std::invalid_argument if the vertex is not in the graph. */
    size_ Degree(V v) const {
      if (!HasVertex(v))
        throw std::invalid_argument("Vertex does not exist.");
      return OutDegree(GetVertexDescriptor(v));
    }
    
    /*! \brief Indegree of a vertex.
     *  \details The indegree of a vertex is the number of edges entering the
     *  vertex. For an undirected graph, this is equivalent to Degree(V) const.
     *  \param v the vertex to get indegree of.
     *  \return the indegree of the vertex.
     *  \throw std::invalid_argument of the vertex is not in the graph. */
    size_ InDegree(V v) const {
      if (!HasVertex(v))
        throw std::invalid_argument("Vertex does not exist.");
      if (D::is_directed)
        return InDegree(GetVertexDescriptor(v));
      return OutDegree(GetVertexDescriptor(v));
    }
    
    /*! \brief Get the neighbouring vertices of a vertex.
     *  \details The neighbours of a vertex are those for which the edge v -> u
     *  exists within the graph.
     *  \param[in] v the vertex to get the neighbours of.
     *  \param[out] nbrs the vector where the list of neighbours will be set.
     *  \return the number of vertices added to the vector.
     *  \throw std::invalid_argument if the vertex is not in the graph. */
    size_ GetNeighbours(V v, std::vector<V>& nbrs) const {
      if (!HasVertex(v))
        throw std::invalid_argument("Vertex does not exist.");
      VertType v_ = GetVertexDescriptor(v);
      nbrs.clear();
      nbrs.reserve(OutDegree(v_));
      NbrsIter begin, end;
      std::tie(begin, end) = boost::adjacent_vertices(v_, *_graph);
      for (; begin != end; ++begin) nbrs.push_back(_verts.right.at(*begin));
      return nbrs.size();
    }
    
    /*! \brief Get the predecessor vertices of a vertex.
     *  \details The predecessors of a vertex are those for which the edge
     *  u -> v exists within the graph. For an undirected graph, this is
     *  equivalent to the neighbours.
     *  \param[in] v the vertex to get the predecessors of.
     *  \param[out] pres the vector where the list of predecessors will be set.
     *  \return the number of vertices added to the vector.
     *  \throw std::invalid_argument if the vertex is not in the graph. */
    size_ GetPredecessors(V v, std::vector<V>& pres) const {
      if (!HasVertex(v))
        throw std::invalid_argument("Vertex does not exist.");
      VertType v_ = GetVertexDescriptor(v);
      pres.clear();
      pres.reserve(OutDegree(v_));
      PredIter begin, end;
      std::tie(begin, end) = boost::inv_adjacent_vertices(v_, *_graph);
      for (; begin != end; ++begin) pres.push_back(_verts.right.at(*begin));
      return pres.size();
    }
    
    /*! \brief Get the two vertices that make up an edge.
     *  \param e the edge to get vertices of.
     *  \return a pair of the two vertices making up the edge.
     *  \throw std::invalid_argument if the edge is not in the graph. */
    std::pair<V, V> GetVertices(E e) const {
      if (!HasEdge(e))
        throw std::invalid_argument("Edge does not exist.");
      EdgeType e_ = GetEdgeDescriptor(e);
      VertType u = boost::source(e_, *_graph);
      VertType v = boost::target(e_, *_graph);
      return std::make_pair(_verts.right.at(u), _verts.right.at(v));
    }
    
    /*! \brief Get the vertices of the graph.
     *  \param[out] verts the vector where the list of vertices will be set.
     *  \return the number of vertices added to the vector. */
    size_ GetVertices(std::vector<V>& verts) const {
      verts.clear();
      verts.reserve(NumVertices());
      VertIter begin, end;
      std::tie(begin, end) = boost::vertices(*_graph);
      for (; begin != end; ++begin) verts.push_back(_verts.right.at(*begin));
      return verts.size();
    }
    
    /*! \brief Get the edges of the graph.
     *  \param[out] edges the vector where the list of edges will be set.
     *  \return the number of edges added to the vector. */
    size_ GetEdges(std::vector<V>& edges) const {
      edges.clear();
      edges.reserve(NumEdges());
      EdgeIter begin, end;
      std::tie(begin, end) = boost::edges(*_graph);
      for (; begin != end; ++begin) edges.push_back(_edges.right.at(*begin));
      return edges.size();
    }
    
    /*! \brief Get the edge between two vertices.
     *  \param u, v vertices to get the edge between.
     *  \return the edge between the two vertces.
     *  \throw std::invalid_argument if either of the vertices is not in the
     *  graph or there is no edge between them. */
    E GetEdge(V u, V v) const {
      if (!HasVertex(u) || !HasVertex(v))
        throw std::invalid_argument("Vertex does not exist.");
      auto e = boost::edge(GetVertexDescriptor(u), GetVertexDescriptor(v),
                           *_graph);
      if (!e.second)
        throw std::invalid_argument("No such edge exists.");
      return _edges.right.at(e.first);
    }
    
    /*! \brief Get the source vertex of an edge.
     *  \param e the edge to get the source of.
     *  \return the source vertex of the edge.
     *  \throw std::invalid_argument if the edge is not in the graph. */
    V GetSource(E e) const {
      if (!HasEdge(e))
        throw std::invalid_argument("Edge does not exist.");
      return _verts.right.at(boost::source(_edges.left.at(e), *_graph));
    }
    
    /*! \brief Get the target vertex of an edge.
     *  \param e the edge to get the target of.
     *  \return the target vertex of the edge.
     *  \throw std::invalid_argument if the edge is not in the graph. */
    V GetTarget(E e) const {
      if (!HasEdge(e))
        throw std::invalid_argument("Edge does not exist.");
      return _verts.right.at(boost::target(_edges.left.at(e), *_graph));
    }
    
  private:
    /*! \brief Get vertex descriptor of a vertex.
     *  \param v vertex to search for.
     *  \return vertex descriptor of the vertex. Behaviour is undefined if
     *  vertex does not exist in graph. */
    VertType GetVertexDescriptor(V v) { return _verts.left.at(v); }
    
    /*! \brief Get edge descriptor of an edge.
     *  \param e edge to search for.
     *  \return pedge descriptor of the edge. Behaviour is undefined if edge
     *  does not exist in graph. */
    EdgeType GetEdgeDescriptor(E e) { return _edges.left.at(e); }
    
    /*! \brief Outdegree of a vertex.
     *  \param v vertex descriptor to get outdegree of.
     *  \return the degree of the given vertex descriptor. */
    size_ OutDegree(VertType v) const { return boost::out_degree(v, *_graph); }
    
    /*! \brief Indegree of a vertex.
     *  \param v vertex descriptor to get indegree of.
     *  \return the indegree of the given vertex descriptor. */
    size_ InDegree(VertType v) const { return boost::in_degree(v, *_graph); }
  };
  
}

#endif /* INDIGOX_GRAPH_BASE_HPP */
