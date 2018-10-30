/*! \file base_graph.hpp */
#ifndef INDIGOX_GRAPH_BASE_HPP
#define INDIGOX_GRAPH_BASE_HPP

#include <cstdint>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "../utils/fwd_declares.hpp"
#include "../utils/modifable_object.hpp"
#include "../utils/simple_bimap.hpp"

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
  
  //! \brief Type for applying a numerical label to a vertex or edge.
  struct GraphLabel{
    //! \brief Union to allow different label types in the same memory space.
    union {
      //! \brief Label used by the connected components algorithm.
      int32_t component;
      //! \brief Label used for (sub-)graph isomorphism.
      uint64_t isomorphism;
    };
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
            class S,
            class D=Undirected,
            class VP=GraphLabel,
            class EP=GraphLabel>
  class BaseGraph : public utils::ModifiableObject {
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    //! \brief Friendship allows graph algorithm access to internals
    friend struct algorithm::access;
  public:
    //! \brief Type of the underlying boost graph.
    using graph_type = boost::adjacency_list<boost::setS,      // Edge container
                                             boost::listS,     // Vertex container
                                 typename D::is_directed_t,    // Directed nature
                                             VP,         // Vertex Properties
                                             EP>;        // Edge Properties
    
    
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
    using NbrsContain = std::map<V, VertContain>;
    //! \brief Type for storing components
    using ComponentContain = std::vector<VertContain>;
    //! \brief Type for storing vertex cycles
    using CycleVertContain = std::vector<VertContain>;
    //! \brief Type for storing edge cycles
    using CycleEdgeContain = std::vector<EdgeContain>;
    using SubgraphType = S;
    
  protected:
    //! \brief Underlying boost graph.
    graph_type _g;
    //! \brief Map vertices to their descriptors.
    VertMap _vm;
    //! \brief Map edges to their descriptors.
    EdgeMap _em;
    //! \brief Container for giving access to all vertices
    VertContain _v;
    //! \brief Container for giving access to all edges
    EdgeContain _e;
    //! \brief Container for predecessors of a vertex (v such that edge u->v exists)
    NbrsContain _pre;
    //! \brief Container for successors of a vertex. Only used in directed graphs
    NbrsContain _suc;
    //! \brief Components container
    ComponentContain _comp_cache;
    //! \brief State when components were last calculated
    utils::ModifiableObject::State _comp_state;
    //! \brief Cyclic vertices container
    VertContain _vcyclic_cache;
    //! \brief Cyclic edges container
    EdgeContain _ecyclic_cache;
    //! \brief Cycles container
    CycleEdgeContain _cycles_cache;
    //! \brief State when cycles were last calculated
    utils::ModifiableObject::State _cycle_state;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t) const;

    template <typename Archive>
    void load(Archive& archive, const uint32_t);
    
  public:
    //! \brief Default constructor
    BaseGraph() : _g(), _comp_state(0), _cycle_state(0) { }
    
  protected:
    // Modification methods protected.
    void Clear();
    
    /*! \brief Add a new vertex to the graph.
     *  \details It is the callers responsability to ensure that the vertex
     *  added is not already part of the graph. If it is, a mismatch between
     *  the vertices in the graph and the what the _verts member thinks are in
     *  the graph may arise.
     *  \param v the vertex to add. */
    void AddVertex(const V& v);
    
    /*! \brief Remove a vertex from the graph.
     *  \details Removing a vertex also removes all edges incident on it. It is
     *  the callers responsibility to ensure that the vertex removed is within
     *  the graph.
     *  \param v the vertex to remove. */
    void RemoveVertex(const V& v);
    
    /*! \brief Add a new edge to the graph.
     *  \details Vertex u is used as the source and vertex v as the target. If
     *  u and/or v are not already part of the graph, they are added. It is the
     *  callers responsibility to ensure that the edge is not part of the graph.
     *  \param u, v vertices the edge is between.
     *  \param e the edge. */
    void AddEdge(const V& u, const V& v, const E& e);
    
    /*! \brief Remove an edge from the graph.
     *  \details It is the callers responsibility to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to remove. */
    void RemoveEdge(const E& e);
    
    /*! \brief Remove an edge from the graph.
     *  \details It is the callers responsibility to ensure that there is an
     *  edge between u and v to remove.
     *  \param u, v vertices to remove an edge from between. */
    void RemoveEdge(const V& u, const V& v);
    
  public:
    /*! \brief Is the vertex in the graph.
     *  \param v vertex to search for.
     *  \return if the requested vertex is contained in the graph or not. */
    bool HasVertex(const V& v) const;
    
    /*! \brief Is the edge in the graph.
     *  \param e edge to search for.
     *  \return if the requested edge is contained in the graph or not. */
    bool HasEdge(const E& e) const;
    
    /*! \brief Does an edge exist between two vertices.
     *  \param u, v vertices to check between.
     *  \return if there is an edge between the two vertices. */
    bool HasEdge(const V& u, const V& v) const;
    
    /*! \brief Number of vertices in the graph.
     *  \return the number of vertices in the graph. */
    int64_t NumVertices() const;
    
    /*! \brief Number of edges in the graph.
     *  \return the number of edges in the graph. */
    int64_t NumEdges() const;
    
    /*! \brief Degree of a vertex.
     *  \details In the case of a directed graph, the degree of a vertex is the
     *  number of edges leaving the vertex.
     *  \param v the vertex to get the degree of.
     *  \return pair of the degree of the vertex and if it is valid. */
    int64_t Degree(const V& v) const;
    
    /*! \brief Indegree of a vertex.
     *  \details The indegree of a vertex is the number of edges entering the
     *  vertex. For an undirected graph, this is equivalent to Degree(V) const.
     *  It is the callers responsibilty to ensure that the vertex is a part of
     *  the graph.
     *  \param v the vertex to get indegree of.
     *  \return pair of the indegree of the vertex and if it is valid. */
    int64_t InDegree(const V& v) const;
    
    /*! \brief Get the neighbouring vertices of a vertex.
     *  \details The neighbours of a vertex are those for which the edge v -> u
     *  exists within the graph. It is the callers responsibilty to ensure that
     *  the vertex is a part of the graph.
     *  \param v the vertex to get the neighbours of.
     *  \param[out] nbrs the vector where the list of neighbours will be set.
     *  The vector is cleared before any neighbouring vertices are added to it.
     *  \return if the vector has been populated or not. */
    const VertContain& GetNeighbours(const V& v);
    
    /*! \brief Get the predecessor vertices of a vertex.
     *  \details The predecessors of a vertex are those for which the edge
     *  u -> v exists within the graph. For an undirected graph, this is
     *  equivalent to the neighbours. It is the callers responsibilty to ensure
     *  that the vertex is a part of the graph.
     *  \param[in] v the vertex to get the predecessors of.
     *  \param[out] pres the vector where the list of predecessors will be set.
     *  The vector is cleared before any predecessing vertices are added to it.
     *  \return if the vector has been populated or not. */
    const VertContain& GetPredecessors(const V& v);
    
    const VertContain& GetSuccessors(const V& v);
    
    /*! \brief Get the two vertices that make up an edge.
     *  \details It is the callers responsibility to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get vertices of.
     *  \return a pair of a pair of the two vertices making up the edge and if
     *  they are valid. */
    std::pair<V, V> GetVertices(const E& e) const;
    
    /*! \brief Get the vertices of the graph.
     *  \param[out] verts the vector where the list of vertices will be set.
     *  The vector is cleared before any vertices are added to it.
     *  \return the number of vertices added to the vector. */
    const VertContain& GetVertices() const;
    
    /*! \brief Get the edges of the graph.
     *  \param[out] edges the vector where the list of edges will be set.
     *  The vector is cleared before any edges are added to it.
     *  \return the number of edges added to the vector. */
    const EdgeContain& GetEdges() const;
    
    /*! \brief Get the edge between two vertices.
     *  \details It is the callers responsibilty to ensure that the vertices
     *  are a part of the graph.
     *  \param u, v vertices to get the edge between.
     *  \return a pair of the edge between the two vertces and if it is valid.*/
    E GetEdge(const V& u, const V& v) const;
    
    /*! \brief Get the source vertex of an edge.
     *  \details It is the callers responsibilty to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get the source of.
     *  \return a pair of the source vertex of the edge and if it is valid. */
    V GetSourceVertex(const E& e) const;
    
    /*! \brief Get the target vertex of an edge.
     *  \details It is the callers responsibilty to ensure that the edge is a
     *  part of the graph.
     *  \param e the edge to get the target of.
     *  \return the target vertex of the edge. */
    V GetTargetVertex(const E& e) const;
    
    /*! \brief Determine if the graph is connected
     *  \return if the graph is connected or not. */
    bool IsConnected();
    
    /*! \brief Get the number of connected components of the graph.
     *  \return the number of cnnected components of the graph. */
    int64_t NumConnectedComponents();
    
    /*! \brief Get the connected components of the graph.
     *  \return reference to the connected components of the graph. */
    const ComponentContain& GetConnectedComponents();
    
    /*! \brief Determine if a vertex of this graph is cyclic
     *  \param v the vertex to check if in a cycle
     *  \return if the vertex is in a cycle or not. */
    bool IsCyclic(const V& v);
    
    /*! \brief Determine if an edge of this graph is cycle.
     *  \param e the edge to check if in a cycle.
     *  \return if the edge is in a cycle or not. */
    bool IsCyclic(const E& e);
    
    /*! \brief Get the cycles of the graph
     *  \return the cycles of the graph. */
    const CycleEdgeContain& GetCycles();
    
    int64_t NumCycles();
    
    virtual S Subgraph(std::vector<V>& verts) = 0;
    virtual S Subgraph(std::vector<V>& verts, std::vector<E>& edges) = 0;
    
  private:
    /*! \brief Get vertex descriptor of a vertex.
     *  \param v vertex to search for.
     *  \return vertex descriptor of the vertex. */
    VertType GetDescriptor(const V& v) const;
    
    /*! \brief Get the vertex of a vertex descriptor.
     *  \param v the vertex descriptor to search for.
     *  \return the vertex associated with the vertex descriptor. */
    V GetV(VertType v) const;
    
    /*! \brief Get edge descriptor of an edge.
     *  \param e edge to search for.
     *  \return edge descriptor of the edge. */
    EdgeType GetDescriptor(const E& e) const;
    
    /*! \brief Get the edge of an edge descriptor.
     *  \param e the edge descriptor to search for.
     *  \return the edge associated with the edge descriptor. */
    E GetE(EdgeType e) const;
    
    /*! \brief Outdegree of a vertex.
     *  \param v vertex descriptor to get outdegree of.
     *  \return the degree of the given vertex descriptor. */
    int64_t OutDegree(VertType v) const;
    
    /*! \brief Indegree of a vertex.
     *  \param v vertex descriptor to get indegree of.
     *  \return the indegree of the given vertex descriptor. */
    int64_t InDegree(VertType v) const;
  };
  
}

#include "base_graph_impl.hpp"

#endif /* INDIGOX_GRAPH_BASE_HPP */
