/*! \file base_graph.hpp */
#ifndef INDIGOX_GRAPH_BASE_HPP
#define INDIGOX_GRAPH_BASE_HPP

#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <indigox/utils/numerics.hpp>
#include <indigox/utils/simple_bimap.hpp>

//! Namespace for all graph related functionality.
namespace indigox::graph {
  
  //! \brief Type for specifying that a graph is directed.
  struct Directed {
    //! \brief Underlying boost type of a directed graph.
    using is_directed_t = boost::bidirectionalS;
    //! \brief Boolean that the type is directed.
    static constexpr bool is_directed = true;
  };
  
  //! \brief Type for specifying that a graph is undirected.
  struct Undirected {
    //! \brief Underlying boost type of an undirected graph.
    using is_directed_t = boost::undirectedS;
    //! \brief Boolean that the type is not directed.
    static constexpr bool is_directed = false;
  };
  
  template <class T>
  class Component : public std::vector<T> {
    using std::vector<T>::vector;
  };
  
  //! \brief Type for applying a numerical label to a vertex or edge.
  struct GraphLabel{
    //! \brief Union to allow different label types in the same memory space.
    union {
      //! \brief Label used by the connected components algorithm.
      int_ component;
      //! \brief Label used for (sub-)graph isomorphism.
      ulong_ isomorphism;
      //! \brief An integer label.
      int_ ilabel;
      //! \brief A floating point label
      float_ flabel;
    };
  };
  
  struct access {
    template <class graph_type, class T>
    inline static graph_type& member_graph(const std::shared_ptr<T>& t) {
      return t._g;
    }
  };
  
  /*! \brief Template base class for all graphs used in the indigoX library.
   *  \tparam V type of the graph vertices.
   *  \tparam E type of the graph edges.
   *  \tparam D type indicating the directed nature of the graph. Defaults to
   *  Undirected.
   *  \tparam VertProp type of the label for a vertex. Defaults to GraphLabel.
   *  \tparam EdgeProp type of the label for an edge. Defaults to GraphLabel.
   *  \details Graph classes should contain an IXGraphBase member, not inherit
   *  from it. Vertex and Edge types are purely interacted with though the use
   *  of raw pointers. The IXGraphBase class will never deallocate any memory
   *  passed to it. */
  template <class V,
            class E,
            class D=Undirected,
            class VertProp=GraphLabel,
            class EdgeProp=GraphLabel>
  class IXGraphBase final {
  private:
    //! \brief Type of the underlying boost graph.
    using graph_t = boost::adjacency_list<boost::setS,      // Edge container
                                          boost::listS,     // Vertex container
                              typename D::is_directed_t,    // Directed nature
                                          VertProp,         // Vertex Properties
                                          EdgeProp>;        // Edge Properties
    
    
    // May need to replace graph_t:: with boost::graph_traits<graph_t>::
    //! \brief Type of the graph vertex descriptor.
    using VertType = typename graph_t::vertex_descriptor;
    //! \brief Type for iterator over graph vertex descriptors.
    using VertIter = typename graph_t::vertex_iterator;
    //! \brief Type for iterator over neighbours of vertex descriptor.
    using NbrsIter = typename graph_t::adjacency_iterator;
    //! \brief Type for iterator over predecessors of a vertex descriptor.
    using PredIter = typename graph_t::inv_adjacency_iterator;
    //! \brief Type of the graph edge descriptor.
    using EdgeType = typename graph_t::edge_descriptor;
    //! \brief Type for iterator over edges.
    using EdgeIter = typename graph_t::edge_iterator;
    
    //! \brief Type for bidirectional mapping of V to vertex descriptor type.
    using VertMap = indigox::utils::SimpleBiMap<V*, VertType>;
    //! \brief Type for bidirectional mapping of E to edge descriptor type.
    using EdgeMap = indigox::utils::SimpleBiMap<E*, EdgeType>;
    //! \brief Friendship allows algorithms access to the underlying boost graph.
    friend struct access;
    
  private:
    //! \brief Underlying boost graph.
    graph_t _g;
    //! \brief Map vertices to their descriptors.
    VertMap _verts;
    //! \brief Map edges to their descriptors.
    EdgeMap _edges;
    
  public:
    //! \brief Default constructor
    IXGraphBase() : _g() { }
    
    /*! \brief Add a new vertex to the graph.
     *  \details It is the callers responsability to ensure that the vertex
     *  added is not already part of the graph. If it is, a mismatch between
     *  the vertices in the graph and the what the _verts member thinks are in
     *  the graph may arise.
     *  \param v the vertex to add. */
    void AddVertex(V* v) {
      VertType v_ = boost::add_vertex(VertProp(), _g);
      _verts.insert(v, v_);
    }
    
    /*! \brief Remove a vertex from the graph.
     *  \details Removing a vertex also removes all edges incident on it. It is
     *  the callers responsibility to ensure that the vertex removed is within
     *  the graph.
     *  \param v the vertex to remove. */
    void RemoveVertex(V* v) {
      VertType v_ = GetDescriptor(v);
      // Remove adjacent edges
      NbrsIter vi, vi_end;
      std::tie(vi, vi_end) = boost::adjacent_vertices(v_, _g);
      for (; vi != vi_end; ++vi) {
        _edges.erase(boost::edge(v_, *vi, _g).first);
      }
      // Remove incident edges of directed graphs
      if (D::is_directed) {
        PredIter vp, vp_end;
        std::tie(vp, vp_end) = boost::inv_adjacent_vertices(v_, _g);
        for (; vp != vp_end; ++vp) {
          _edges.erase(boost::edge(*vp, v_, _g).first);
        }
      }
      // Remove the vertex
      _verts.erase(v_);
      boost::clear_vertex(v_, _g);
      boost::remove_vertex(v_, _g);
    }
    
    /*! \brief Add a new edge to the graph.
     *  \details Vertex u is used as the source and vertex v as the target. If
     *  u and/or v are not already part of the graph, they are added. It is the
     *  callers responsibility to ensure that the edge is not part of the graph.
     *  \param u, v vertices the edge is between.
     *  \param e the edge. */
    void AddEdge(V* u, V* v, E* e) {
      if (!HasVertex(u)) AddVertex(u);
      if (!HasVertex(v)) AddVertex(v);
      VertType u_ = GetDescriptor(u);
      VertType v_ = GetDescriptor(v);
      
      EdgeType e_ = boost::add_edge(u_, v_, EdgeProp(), _g).first;
      _edges.insert(e, e_);
    }
    
    /*! \brief Remove an edge from the graph.
     *  \details It is the callers responsibility to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to remove. */
    void RemoveEdge(E* e) {
      EdgeType e_ = GetDescriptor(e);
      _edges.erase(e_);
      boost::remove_edge(e_, _g);
    }
    
    /*! \brief Remove an edge from the graph.
     *  \details It is the callers responsibility to ensure that there is an
     *  edge between u and v to remove.
     *  \param u, v vertices to remove an edge from between. */
    void RemoveEdge(V* u, V* v) {
      VertType u_ = GetDescriptor(u);
      VertType v_ = GetDescriptor(v);
      EdgeType e = boost::edge(u_, v_, _g).first;
      _edges.erase(e);
      boost::remove_edge(e, _g);
    }
    
    /*! \brief Is the vertex in the graph.
     *  \param v vertex to search for.
     *  \return if the requested vertex is contained in the graph or not. */
    bool HasVertex(V* v) const {
      return _verts.left.find(v) != _verts.left.end();
    }
    
    /*! \brief Is the edge in the graph.
     *  \param e edge to search for.
     *  \return if the requested edge is contained in the graph or not. */
    bool HasEdge(E* e) const {
      return _edges.left.find(e) != _edges.left.end();
    }
    
    /*! \brief Does an edge exist between two vertices.
     *  \param u, v vertices to check between.
     *  \return if there is an edge between the two vertices. */
    bool HasEdge(V* u, V* v) const {
      if (!HasVertex(u) || !HasVertex(v)) return false;
      VertType u_ = GetDescriptor(u);
      VertType v_ = GetDescriptor(v);
      return boost::edge(u_, v_, _g).second;
    }
    
    /*! \brief Number of vertices in the graph.
     *  \return the number of vertices in the graph. */
    size_ NumVertices() const { return boost::num_vertices(_g); }
    
    /*! \brief Number of edges in the graph.
     *  \return the number of edges in the graph. */
    size_ NumEdges() const { return boost::num_edges(_g); }
    
    //! \brief Removes all edges and vertices from the graph.
    void Clear() {
      _g.clear();
      _verts.clear();
      _edges.clear();
    }
    
    /*! \brief Degree of a vertex.
     *  \details In the case of a directed graph, the degree of a vertex is the
     *  number of edges leaving the vertex. It is the callers responsibilty to
     *  ensure that the vertex is a prt of the graph.
     *  \param v the vertex to get the degree of.
     *  \return pair of the degree of the vertex and if it is valid. */
    size_ Degree(V* v) const { return OutDegree(GetDescriptor(v)); }
    
    /*! \brief Indegree of a vertex.
     *  \details The indegree of a vertex is the number of edges entering the
     *  vertex. For an undirected graph, this is equivalent to Degree(V) const.
     *  It is the callers responsibilty to ensure that the vertex is a part of
     *  the graph.
     *  \param v the vertex to get indegree of.
     *  \return pair of the indegree of the vertex and if it is valid. */
    size_ InDegree(V* v) const {
      if (D::is_directed) return InDegree(GetDescriptor(v));
      return OutDegree(GetDescriptor(v));
    }
    
    /*! \brief Get the neighbouring vertices of a vertex.
     *  \details The neighbours of a vertex are those for which the edge v -> u
     *  exists within the graph. It is the callers responsibilty to ensure that
     *  the vertex is a part of the graph.
     *  \param v the vertex to get the neighbours of.
     *  \param[out] nbrs the vector where the list of neighbours will be set.
     *  The vector is cleared before any neighbouring vertices are added to it.
     *  \return if the vector has been populated or not. */
    void GetNeighbours(V* v, std::vector<V*>& nbrs) const {
      VertType v_ = GetDescriptor(v);
      nbrs.clear(); nbrs.reserve(OutDegree(v_));
      NbrsIter begin, end;
      std::tie(begin, end) = boost::adjacent_vertices(v_, _g);
      for (; begin != end; ++begin) nbrs.push_back(_verts.right.at(*begin));
    }
    
    /*! \brief Get the predecessor vertices of a vertex.
     *  \details The predecessors of a vertex are those for which the edge
     *  u -> v exists within the graph. For an undirected graph, this is
     *  equivalent to the neighbours. It is the callers responsibilty to ensure
     *  that the vertex is a part of the graph.
     *  \param[in] v the vertex to get the predecessors of.
     *  \param[out] pres the vector where the list of predecessors will be set.
     *  The vector is cleared before any predecessing vertices are added to it.
     *  \return if the vector has been populated or not. */
    void GetPredecessors(V* v, std::vector<V*>& pres) const {
      if (!D::is_directed) {
        GetNeighbours(v, pres);
        return;
      }
      VertType v_ = GetDescriptor(v);
      pres.clear(); pres.reserve(InDegree(v_));
      PredIter begin, end;
      std::tie(begin, end) = boost::inv_adjacent_vertices(v_, _g);
      for (; begin != end; ++begin) pres.push_back(_verts.right.at(*begin));
    }
    
    /*! \brief Get the two vertices that make up an edge.
     *  \details It is the callers responsibility to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get vertices of.
     *  \return a pair of a pair of the two vertices making up the edge and if
     *  they are valid. */
    std::pair<V*, V*> GetVertices(E* e) const {
      EdgeType e_ = GetDescriptor(e);
      VertType u = boost::source(e_, _g);
      VertType v = boost::target(e_, _g);
      return {_verts.right.at(u), _verts.right.at(v)};
    }
    
    /*! \brief Get the vertices of the graph.
     *  \param[out] verts the vector where the list of vertices will be set.
     *  The vector is cleared before any vertices are added to it.
     *  \return the number of vertices added to the vector. */
    size_ GetVertices(std::vector<V*>& verts) const {
      verts.clear(); verts.reserve(NumVertices());
      auto begin = _verts.left.begin();
      auto end = _verts.left.end();
      for (; begin != end; ++begin) verts.push_back(begin->first);
      return verts.size();
    }
    
    /*! \brief Get the edges of the graph.
     *  \param[out] edges the vector where the list of edges will be set.
     *  The vector is cleared before any edges are added to it.
     *  \return the number of edges added to the vector. */
    size_ GetEdges(std::vector<E*>& edges) const {
      edges.clear(); edges.reserve(NumEdges());
      EdgeIter begin, end;
      std::tie(begin, end) = boost::edges(_g);
      for (; begin != end; ++begin) edges.push_back(_edges.right.at(*begin));
      return edges.size();
    }
    
    /*! \brief Get the edge between two vertices.
     *  \details It is the callers responsibilty to ensure that the vertices
     *  are a part of the graph.
     *  \param u, v vertices to get the edge between.
     *  \return a pair of the edge between the two vertces and if it is valid.*/
    E* GetEdge(V* u, V* v) const {
      auto e = boost::edge(GetDescriptor(u), GetDescriptor(v), _g);
      return _edges.right.at(e.first);
    }
    
    /*! \brief Get the source vertex of an edge.
     *  \details It is the callers responsibilty to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get the source of.
     *  \return a pair of the source vertex of the edge and if it is valid. */
    V* GetSource(E* e) const {
      return _verts.right.at(boost::source(_edges.left.at(e), _g));
    }
    
    /*! \brief Get the target vertex of an edge.
     *  \details It is the callers responsibilty to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get the target of.
     *  \return the target vertex of the edge. */
    V* GetTarget(E* e) const {
      return _verts.right.at(boost::target(_edges.left.at(e), _g));
    }
    
    /*! \brief Get the connected components of the graph.
     *  \param[out] components where the connected components will be written to.
     *  \return the number of connected components. */
//    size_ ConnectedComponents(std::vector<std::vector<V*>>& components) {
//      static_assert(!D::is_directed, "Requires an undirected graph.");
//      size_ num = __connected_components_worker();
//      components.clear();
//      components.assign(num, std::vector<V*>());
//      for (auto v : _verts.right)
//        components[(*_graph)[v.first].ilabel].push_back(v.second);
//      return num;
//    }
    
    /*! \brief Number of connected components of the graph.
     *  \return the number of connected components. */
//    size_ NumConnectedComponents() {
//      static_assert(!D::is_directed, "Requires an undirected graph.");
//      return __connected_components_worker();
//    }
    
  private:
    /*! \brief Get vertex descriptor of a vertex.
     *  \param v vertex to search for.
     *  \return vertex descriptor of the vertex. */
    VertType GetDescriptor(V* v) const { return _verts.left.at(v); }
    
    /*! \brief Get the vertex of a vertex descriptor.
     *  \param v the vertex descriptor to search for.
     *  \return the vertex associated with the vertex descriptor. */
    V* GetVertex(VertType v) const { return _verts.right.at(v); }
    
    /*! \brief Get edge descriptor of an edge.
     *  \param e edge to search for.
     *  \return edge descriptor of the edge. */
    EdgeType GetDescriptor(E* e) const { return _edges.left.at(e); }
    
    /*! \brief Get the edge of an edge descriptor.
     *  \param e the edge descriptor to search for.
     *  \return the edge associated with the edge descriptor. */
    E* GetEdge(EdgeType e) const { return _edges.right.at(e); }
    
    /*! \brief Outdegree of a vertex.
     *  \param v vertex descriptor to get outdegree of.
     *  \return the degree of the given vertex descriptor. */
    size_ OutDegree(VertType v) const { return boost::out_degree(v, _g); }
    
    /*! \brief Indegree of a vertex.
     *  \param v vertex descriptor to get indegree of.
     *  \return the indegree of the given vertex descriptor. */
    size_ InDegree(VertType v) const { return boost::in_degree(v, _g); }
    
    //! \cond
    /*! \brief Worker method for calculating connected components.
     *  \details Each vertex is labeled with the ID of the component it is
     *  a part of. */
//    size_ __connected_components_worker() {
//      using namespace boost;
//      using VertIdxMap = std::map<VertType, size_>;
//      VertIdxMap data;
//      associative_property_map<VertIdxMap> indexMap(data);
//      auto verts = vertices(*_graph);
//      for (size_ i = 0; verts.first != verts.second; ++verts.first, ++i)
//        put(indexMap, *verts.first, i);
//
//      auto num = connected_components(*_graph,
//                                      get(&VertProp::component, *_graph),
//                                      vertex_index_map(indexMap));
//      return static_cast<size_>(num);
//    }
    //! \endcond
  };
  
}

#endif /* INDIGOX_GRAPH_BASE_HPP */
