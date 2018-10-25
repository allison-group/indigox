/*! \file base_graph.hpp */
#ifndef INDIGOX_GRAPH_BASE_HPP
#define INDIGOX_GRAPH_BASE_HPP

#include <cstdint>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "../utils/simple_bimap.hpp"
#include "../utils/triple.hpp"

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
      int32_t component;
      //! \brief Label used for (sub-)graph isomorphism.
      uint64_t isomorphism;
      //! \brief An integer label.
      int32_t ilabel;
      //! \brief A floating point label
      double flabel;
      //! \brief A colour for graph algorithms
      boost::default_color_type colour;
    };
  };
  
  struct access {
    template <class T>
    inline static typename T::graph_type& member_graph(const std::shared_ptr<T>& t) {
      return t->_g;
    }
    
    template <class T>
    inline static typename T::graph_type& member_graph(T& t) {
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
  class BaseGraph {
    //! \brief Friendship allows serialisation
    friend class cereal::access;
  public:
    //! \brief Type of the underlying boost graph.
    using graph_type = boost::adjacency_list<boost::setS,      // Edge container
                                             boost::listS,     // Vertex container
                                 typename D::is_directed_t,    // Directed nature
                                             VertProp,         // Vertex Properties
                                             EdgeProp>;        // Edge Properties
    
    
    // May need to replace graph_t:: with boost::graph_traits<graph_t>::
    //! \brief Type of the graph vertex descriptor.
    using VertType = typename graph_type::vertex_descriptor;
    //! \brief Type for iterator over graph vertex descriptors.
    using VertIter = typename graph_type::vertex_iterator;
    //! \brief Type for iterator over neighbours of vertex descriptor.
    using NbrsIter = typename graph_type::adjacency_iterator;
    //! \brief Type for iterator over predecessors of a vertex descriptor.
    using PredIter = typename graph_type::inv_adjacency_iterator;
    //! \brief Type of the graph edge descriptor.
    using EdgeType = typename graph_type::edge_descriptor;
    //! \brief Type for iterator over edges.
    using EdgeIter = typename graph_type::edge_iterator;
    
    //! \brief Type for bidirectional mapping of V to vertex descriptor type.
    using VertMap = indigox::utils::SimpleBiMap<V, VertType>;
    //! \brief Type for bidirectional mapping of E to edge descriptor type.
    using EdgeMap = indigox::utils::SimpleBiMap<E, EdgeType>;
    //! \brief Type for storing vertices.
    using VertContain = std::vector<V>;
    //! \brief Type for storing edges
    using EdgeContain = std::vector<E>;
    //! \brief Type for storing neighbours
    using NbrsContain = std::map<V, std::vector<V>>;
    //! \brief Friendship allows algorithms access to the underlying boost graph.
    friend struct access;
    
  public:
    //! \brief Underlying boost graph.
    graph_type _g;
    //! \brief Map vertices to their descriptors.
    VertMap _vm;
    //! \brief Map edges to their descriptors.
    EdgeMap _em;
    //! \brief Container for giving iterator access to all vertices
    VertContain _v;
    //! \brief Container for giving iterator access to all edges
    EdgeContain _e;
    //! \brief Container for predecessors of a vertex (v such that edge u->v exists)
    NbrsContain _pre;
    //! \brief Container for successors of a vertex. Only used in directed graphs
    NbrsContain _suc;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t) const {
      std::vector<stdx::triple<V, V, E>> edges;
      edges.reserve(NumEdges());
      for (const E& e : _e)
        edges.emplace_back(GetSourceVertex(e), GetTargetVertex(e), e);

      archive(INDIGOX_SERIAL_NVP("vertices", _v),
              INDIGOX_SERIAL_NVP("edges", edges));
    }

    template <typename Archive>
    void load(Archive& archive, const uint32_t) {
      VertContain vertices;
      std::vector<stdx::triple<V, V, E>> edges;
      archive(INDIGOX_SERIAL_NVP("vertices", vertices),
              INDIGOX_SERIAL_NVP("edges", edges));

      // Build the graph
      for (V& v : vertices) AddVertex(v);
      for (auto& e : edges) AddEdge(e.first, e.second, e.third);
    }
    
  public:
    //! \brief Default constructor
    BaseGraph() : _g() { }
    
  protected:
    // Modification methods protected.
    void Clear() {
      _vm.clear();
      _em.clear();
      _v.clear();
      _e.clear();
      _pre.clear();
      _suc.clear();
      _g.clear();
    }
    
    /*! \brief Add a new vertex to the graph.
     *  \details It is the callers responsability to ensure that the vertex
     *  added is not already part of the graph. If it is, a mismatch between
     *  the vertices in the graph and the what the _verts member thinks are in
     *  the graph may arise.
     *  \param v the vertex to add. */
    void AddVertex(const V& v) {
      VertType vboost = boost::add_vertex(VertProp(), _g);
      _vm.insert(v, vboost);
      _v.emplace_back(v);
      _pre.emplace(v, VertContain());
      if (D::is_directed) _suc.emplace(v, VertContain());
    }
    
    /*! \brief Remove a vertex from the graph.
     *  \details Removing a vertex also removes all edges incident on it. It is
     *  the callers responsibility to ensure that the vertex removed is within
     *  the graph.
     *  \param v the vertex to remove. */
    void RemoveVertex(const V& v) {
      VertType vboost = GetDescriptor(v);
      // Remove adjacent edges
      PredIter vi, vi_end;
      std::tie(vi, vi_end) = boost::inv_adjacent_vertices(vboost, _g);
      for (; vi != vi_end; ++vi) {
        V u = GetV(*vi);
        E e = GetE(boost::edge(vboost, *vi, _g).first);
        _em.erase(e);
        _e.erase(std::find(_e.begin(), _e.end(), e));
      }
      // Remove incident edges of directed graphs
      if (D::is_directed) {
        NbrsIter vp, vp_end;
        std::tie(vp, vp_end) = boost::adjacent_vertices(vboost, _g);
        for (; vp != vp_end; ++vp) {
          V u = GetV(*vp);
          E e = GetE(boost::edge(*vp, vboost, _g).first);
          _em.erase(e);
          _e.erase(std::find(_e.begin(), _e.end(), e));
        }
      }
      // Remove the vertex
      _vm.erase(v);
      _v.erase(std::find(_v.begin(), _v.end(), v));
      _pre.erase(v);
      if (D::is_directed) _suc.erase(v);
      boost::clear_vertex(vboost, _g);
      boost::remove_vertex(vboost, _g);
    }
    
    /*! \brief Add a new edge to the graph.
     *  \details Vertex u is used as the source and vertex v as the target. If
     *  u and/or v are not already part of the graph, they are added. It is the
     *  callers responsibility to ensure that the edge is not part of the graph.
     *  \param u, v vertices the edge is between.
     *  \param e the edge. */
    void AddEdge(const V& u, const V& v, const E& e) {
      if (!HasVertex(u)) AddVertex(u);
      if (!HasVertex(v)) AddVertex(v);
      VertType uboost = GetDescriptor(u);
      VertType vboost = GetDescriptor(v);
      EdgeType eboost = boost::add_edge(uboost, vboost, EdgeProp(), _g).first;
      _em.insert(e, eboost);
      _e.emplace_back(e);
    }
    
    /*! \brief Remove an edge from the graph.
     *  \details It is the callers responsibility to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to remove. */
    void RemoveEdge(const E& e) {
      EdgeType eboost = GetDescriptor(e);
      _em.erase(eboost);
      _e.erase(std::find(_e.begin(), _e.end(), e));
      boost::remove_edge(eboost, _g);
    }
    
    /*! \brief Remove an edge from the graph.
     *  \details It is the callers responsibility to ensure that there is an
     *  edge between u and v to remove.
     *  \param u, v vertices to remove an edge from between. */
    void RemoveEdge(const V& u, const V& v) {
      E e = GetEdge(u, v);
      RemoveEdge(e);
    }
    
  public:
    /*! \brief Is the vertex in the graph.
     *  \param v vertex to search for.
     *  \return if the requested vertex is contained in the graph or not. */
    bool HasVertex(const V& v) const {
      return _vm.left.find(v) != _vm.left.end();
    }
    
    /*! \brief Is the edge in the graph.
     *  \param e edge to search for.
     *  \return if the requested edge is contained in the graph or not. */
    bool HasEdge(const E& e) const {
      return _em.left.find(e) != _em.left.end();
    }
    
    /*! \brief Does an edge exist between two vertices.
     *  \param u, v vertices to check between.
     *  \return if there is an edge between the two vertices. */
    bool HasEdge(const V& u, const V& v) const {
      if (!HasVertex(u) || !HasVertex(v)) return false;
      VertType u_ = GetDescriptor(u);
      VertType v_ = GetDescriptor(v);
      return boost::edge(u_, v_, _g).second;
    }
    
    /*! \brief Number of vertices in the graph.
     *  \return the number of vertices in the graph. */
    int64_t NumVertices() const { return boost::num_vertices(_g); }
    
    /*! \brief Number of edges in the graph.
     *  \return the number of edges in the graph. */
    int64_t NumEdges() const { return boost::num_edges(_g); }
    
    /*! \brief Degree of a vertex.
     *  \details In the case of a directed graph, the degree of a vertex is the
     *  number of edges leaving the vertex.
     *  \param v the vertex to get the degree of.
     *  \return pair of the degree of the vertex and if it is valid. */
    int64_t Degree(const V& v) const {
      return HasVertex(v) ? OutDegree(GetDescriptor(v)) : -1;
    }
    
    /*! \brief Indegree of a vertex.
     *  \details The indegree of a vertex is the number of edges entering the
     *  vertex. For an undirected graph, this is equivalent to Degree(V) const.
     *  It is the callers responsibilty to ensure that the vertex is a part of
     *  the graph.
     *  \param v the vertex to get indegree of.
     *  \return pair of the indegree of the vertex and if it is valid. */
    int64_t InDegree(const V& v) const {
      if (!HasVertex(v)) return -1;
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
    const VertContain& GetNeighbours(const V& v) {
      static_assert(!D::is_directed, "Requires an undirected graph.");
      VertType vboost = GetDescriptor(v);
      auto adjis = boost::adjacent_vertices(vboost, _g);
      _pre.at(v).clear();
      for (; adjis.first != adjis.second; ++adjis.first)
        _pre.at(v).emplace_back(GetV(*adjis.first));
      return _pre.at(v);
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
    const VertContain& GetPredecessors(const V& v) {
      static_assert(D::is_directed, "Requires a directed graph.");
      VertType vboost = GetDescriptor(v);
      auto adjis = boost::inv_adjacent_vertices(vboost, _g);
      _pre.at(v).clear();
      for (; adjis.first != adjis.second; ++adjis.first)
        _pre.at(v).emplace_back(GetV(*adjis.first));
      return _pre.at(v);
    }
    
    const VertContain& GetSuccessors(const V& v) {
      static_assert(D::is_directed, "Requires a directed graph.");
      VertType vboost = GetDescriptor(v);
      auto adjis = boost::adjacent_vertices(vboost, _g);
      _suc.at(v).clear();
      for (; adjis.first != adjis.second; ++adjis.first)
        _suc.at(v).emplace_back(GetV(*adjis.first));
      return _suc.at(v);
    }
    
    /*! \brief Get the two vertices that make up an edge.
     *  \details It is the callers responsibility to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get vertices of.
     *  \return a pair of a pair of the two vertices making up the edge and if
     *  they are valid. */
    std::pair<V, V> GetVertices(const E& e) const {
      EdgeType eboost = GetDescriptor(e);
      VertType u = boost::source(eboost, _g);
      VertType v = boost::target(eboost, _g);
      return {GetV(u), GetV(v)};
    }
    
    /*! \brief Get the vertices of the graph.
     *  \param[out] verts the vector where the list of vertices will be set.
     *  The vector is cleared before any vertices are added to it.
     *  \return the number of vertices added to the vector. */
    const VertContain& GetVertices() const { return _v; }
    
    /*! \brief Get the edges of the graph.
     *  \param[out] edges the vector where the list of edges will be set.
     *  The vector is cleared before any edges are added to it.
     *  \return the number of edges added to the vector. */
    const EdgeContain& GetEdges() const { return _e; }
    
    /*! \brief Get the edge between two vertices.
     *  \details It is the callers responsibilty to ensure that the vertices
     *  are a part of the graph.
     *  \param u, v vertices to get the edge between.
     *  \return a pair of the edge between the two vertces and if it is valid.*/
    E GetEdge(const V& u, const V& v) const {
      VertType uboost = GetDescriptor(u);
      VertType vboost = GetDescriptor(v);
      EdgeType eboost = boost::edge(uboost, vboost, _g).first;
      return GetE(eboost);
    }
    
    /*! \brief Get the source vertex of an edge.
     *  \details It is the callers responsibilty to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get the source of.
     *  \return a pair of the source vertex of the edge and if it is valid. */
    V GetSourceVertex(const E& e) const {
      EdgeType eboost = GetDescriptor(e);
      VertType source = boost::source(eboost, _g);
      return GetV(source);
    }
    
    /*! \brief Get the target vertex of an edge.
     *  \details It is the callers responsibilty to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get the target of.
     *  \return the target vertex of the edge. */
    V GetTargetVertex(const E& e) const {
      EdgeType eboost = GetDescriptor(e);
      VertType target = boost::target(eboost, _g);
      return GetV(target);
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
    VertType GetDescriptor(const V& v) const { return _vm.left.at(v); }
    
    /*! \brief Get the vertex of a vertex descriptor.
     *  \param v the vertex descriptor to search for.
     *  \return the vertex associated with the vertex descriptor. */
    V GetV(VertType v) const { return _vm.right.at(v); }
    
    /*! \brief Get edge descriptor of an edge.
     *  \param e edge to search for.
     *  \return edge descriptor of the edge. */
    EdgeType GetDescriptor(const E& e) const { return _em.left.at(e); }
    
    /*! \brief Get the edge of an edge descriptor.
     *  \param e the edge descriptor to search for.
     *  \return the edge associated with the edge descriptor. */
    E GetE(EdgeType e) const { return _em.right.at(e); }
    
    /*! \brief Outdegree of a vertex.
     *  \param v vertex descriptor to get outdegree of.
     *  \return the degree of the given vertex descriptor. */
    size_t OutDegree(VertType v) const { return boost::out_degree(v, _g); }
    
    /*! \brief Indegree of a vertex.
     *  \param v vertex descriptor to get indegree of.
     *  \return the indegree of the given vertex descriptor. */
    size_t InDegree(VertType v) const { return boost::in_degree(v, _g); }
    
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
