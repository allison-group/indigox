/*! \file condensed.hpp */
#ifndef INDIGOX_GRAPH_CONDENSED_HPP
#define INDIGOX_GRAPH_CONDENSED_HPP

#include <memory>
#include <vector>

#include <EASTL/bitset.h>
#include <EASTL/vector_map.h>
#include <EASTL/vector_set.h>

#include "base_graph.hpp"
#include "../utils/common.hpp"
#include "../utils/numerics.hpp"


// Forward declares
namespace indigox::test {
  struct TestCondensedMolecularGraph;
  struct TestCondensedVertex;
  struct TestCondensedEdge;
}
namespace indigox::graph {
  class IXMolecularGraph;
  class IXMGVertex;
  class IXMGEdge;
  
  using MolecularGraph = std::shared_ptr<IXMolecularGraph>;
  using MGVertex = std::shared_ptr<IXMGVertex>;
  using MGEdge = std::shared_ptr<IXMGEdge>;
  
  using _MolecularGraph = std::weak_ptr<IXMolecularGraph>;
  using _MGVertex = std::weak_ptr<IXMGVertex>;
  using _MGEdge = std::weak_ptr<IXMGEdge>;
  
  class IXCondensedMolecularGraph;
  class IXCMGVertex;
  class IXCMGEdge;
  
  //! \brief shared_ptr for normal use of the IXCondensedMolecularGraph class.
  using CondensedMolecularGraph = std::shared_ptr<IXCondensedMolecularGraph>;
  //! \brief shared_ptr for normal use of the IXCMGVertex class.
  using CMGVertex = std::shared_ptr<IXCMGVertex>;
  //! \brief shared_ptr for normal use of the IXCMGEdge class.
  using CMGEdge = std::shared_ptr<IXCMGEdge>;
  
  /*! \brief weak_ptr for non-ownership reference to the IXCondensedMolecularGraph class.
   *  \details Intended for internal use only. */
  using _CondensedMolecularGraph = std::weak_ptr<IXCondensedMolecularGraph>;
  /*! \brief weak_ptr for non-ownership reference to the IXCMGVertex class.
   *  \details Intended for internal use only. */
  using _CMGVertex = std::weak_ptr<IXCMGVertex>;
  /*! \brief weak_ptr for non-ownership reference to the IXCMGEdge class.
   *  \details Intended for internal use only. */
  using _CMGEdge = std::weak_ptr<IXCMGEdge>;
  
  //! \brief type used to store the isomorphism testing mask for IXCMGVertex.
  using VertexIsoMask = eastl::bitset<32, uint64_t>;
  //! \brief type used to store the isomorphism testing mask for IXCMGEdge.
  using EdgeIsoMask = eastl::bitset<16, uint32_t>;
}


namespace indigox::graph {
  
  inline bool operator<(_CMGVertex lhs, _CMGVertex rhs) {
    return lhs.lock() < rhs.lock();
  }
  
  /*! \brief Class for the vertices of an IXCondensedMolecularGraph. */
  class IXCMGVertex : public std::enable_shared_from_this<IXCMGVertex> {
  public:
    //! \brief Enum for the contracted symmetry groups
    enum class ContractedSymmetry {
      Hydrogen, //!< Group for contracted hydrogen atoms.
      Flourine, //!< Group for contracted flourine atoms.
      Chlorine, //!< Group for contracted chlorine atoms.
      Bromine,  //!< Group for contracted bromine atoms.
      Iodine    //!< Group for contracted iodine atoms.
    };
    
  private:
    //! \brief Friendship allows IXCondensedMolecularGraph to add vertices.
    friend class IXCondensedMolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;
    //! \brief Friendship allows for testing
    friend struct indigox::test::TestCondensedVertex;
    
    //! \brief Type used to store condensed vertices
    using CondensedVertex = std::pair<ContractedSymmetry, MGVertex>;
    
  public:
    IXCMGVertex() = delete; // no default constructor
    
    /*! \brief Get the MGVertex associated with this vertex.
     *  \return the associated MGVertex, if it is still alive. */
    inline MGVertex GetSource() const { return _source; }
    
    /*! \brief Get the graph this vertex is part of.
     *  \return the owning graph. */
    inline CondensedMolecularGraph GetGraph() const { return _graph.lock(); }
    
    /*! \brief Get the number of contracted vertices.
     *  \return the number of contracted vertices. */
    inline size_ NumContracted() const { return _contracted.size(); }
    
    /*! \brief Get the number of contracted vertices in the symmetry group.
     *  \return the number of vertices of the given symmetry group contracted. */
    size_ NumContracted(ContractedSymmetry sym) const;
    
    /*! \brief Get the isomorphism testing mask.
     *  \return the isomorphism testing mask. */
    inline VertexIsoMask GetIsomorphismMask() const { return _iso_mask; }
    
    /*! \brief Checks if a given MGVertex is contracted into this vertex.
     *  \param v the vertex to check for.
     *   return if the provided vertex is contracted into this one. */
    bool IsContractedHere(const MGVertex& v) const;
    
    // Not worrying about getting the contracted vertices just yet
    
  private:
    /*! \brief Construct a vertex from an MGVertex.
     *  \param v the MGVertex to associate with this vertex.
     *  \param g the condensed graph this vertex will be a part of. */
    IXCMGVertex(const MGVertex& v, const CondensedMolecularGraph& g);
    
  private:
    //! \brief Source vertex
    MGVertex _source;
    //! \brief Parent graph
    _CondensedMolecularGraph _graph;
    //! \brief Contracted vertices
    eastl::vector_set<CondensedVertex> _contracted;
    //! \brief Isomorphism testing mask
    VertexIsoMask _iso_mask;
  };
  
  /*! \brief Class for the edges of an IXCondensedMolecularGraph. */
  class IXCMGEdge : public std::enable_shared_from_this<IXCMGEdge> {
    //! \brief Friendship allows IXCondensedMolecularGraph to add edges.
    friend class IXCondensedMolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;
    //! \brief Friendship allows for testing
    friend struct indigox::test::TestCondensedEdge;
    
  public:
    IXCMGEdge() = delete; // no default constructor
    
    /*! \brief Get the MGEdge associated with this vertex.
     *  \return the associated MGEdge, if it is atill alive. */
    inline MGEdge GetSource() const { return _source.lock(); }
    
    /*! \brief Get the graph this edge is part of.
     *  \return the owning graph. */
    inline CondensedMolecularGraph GetGraph() const { return _graph.lock(); }
    
    /*! \brief Get the isomorphism testing mask.
     *  \return the isomorphism testing mask. */
    inline EdgeIsoMask GetIsomorphismMask() const { return _iso_mask; }
    
  private:
    /*! \brief Construct an edge from an MGEdge.
     *  \param e the MGEdge to associate with this edge.
     *  \param g the condensed graph this edge will be a part of. */
    IXCMGEdge(const MGEdge& e, const CondensedMolecularGraph& g);
    
  private:
    //! \brief Source edge
    _MGEdge _source;
    //! \brief Parent graph
    _CondensedMolecularGraph _graph;
    //! \brief Isomorhism testing mask
    EdgeIsoMask _iso_mask;
  };
  
  class IXCondensedMolecularGraph:
  public std::enable_shared_from_this<IXCondensedMolecularGraph> {
    //! \brief Friendship allows IXCondensedMolecularGraph to be tested.
    friend struct indigox::test::TestCondensedMolecularGraph;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    //! \brief Friendship allows graph algorithm access to underlying graph
    friend struct indigox::graph::access;
    
    //! \brief Type of the underlying IXGraphBase
    using graph_type = IXGraphBase<IXCMGVertex, IXCMGEdge>;
    //! \brief Container for vertices
    using VertContain = std::vector<CMGVertex>;
    //! \brief Container for edges
    using EdgeContain = std::vector<CMGEdge>;
    //! \brief Container for neighbours of vertices
    using NbrsContain = std::map<CMGVertex, VertContain>;
    //! \brief Container for mapping MGVertex to vertices
    using VertMap = eastl::vector_map<MGVertex, CMGVertex>;
    //! \brief Container for mapping MGEdge to edges
    using EdgeMap = eastl::vector_map<MGEdge, CMGEdge>;
    
  public:
    //! \brief Type of the iterator returned by GetEdges() method.
    using EdgeIter = EdgeContain::const_iterator;
    //! \brief Type of the iterator returned by GetVertices() method.
    using VertIter = VertContain::const_iterator;
    //! \brief Type of the iterator returned by GetNeighbours() method.
    using NbrsIter = NbrsContain::mapped_type::const_iterator;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                      cereal::construct<IXCondensedMolecularGraph>& construct,
                                   const uint32_t version);
    
    // Modification methods are private so that the CMG is a snapshot of the MG
    // at time of creation
    /*! \brief Add an edge to the graph.
     *  \details Adds an edge between previously added vertices. Assumes that
     *  the source and target vertex of the provided MGEdge have already been
     *  added to the graph.
     *  \param e the source MGEdge.
     *  \return shared_ptr tp the newly added edge. */
    CMGEdge AddEdge(const MGEdge& e);
    
    /*! \brief Add a vertex to the graph.
     *  \details It is assumed that the provided source vertex is not viable
     *  for condensing.
     *  \param v the source MGVertex.
     *  \return shared_ptr to the newly added vertex. */
    CMGVertex AddVertex(const MGVertex& v);
    
    // No need for removal methods as snapshot will never need them
    
  public:
    IXCondensedMolecularGraph() = default;  // no default constructor
    
    /*! \brief Construct a condensed molecular graph from a MolecularGraph.
     *  \param g the molecular graph to construct from. */
    IXCondensedMolecularGraph(const MolecularGraph& g);
    
    /*! \brief Get the source MolecularGraph.
     *  \return the molecular graph used to construt this. */
    inline MolecularGraph GetSource() const { return _source; }
    
    /*! \brief The degree of a vertex.
     *  \details If the vertex is not part of the graph, the returned values is
     *  std::numeric_limits<size_>::max().
     *  \param v the vertex to obtain the degree of.
     *  \return the degree of the vertex. */
    size_ Degree(const CMGVertex v) const;
    
    /*! \brief Get the edge between two vertices.
     *  \details If there is no edge between the vertices, or at least one of
     *  the vertices is not part of the graph, the returned edge is null.
     *  \param u, v the vertices to get the edge between.
     *  \return the edge between the two vertices. */
    CMGEdge GetEdge(const CMGVertex u, const CMGVertex v) const;
    
    /*! \brief Get the edge associated with an MGEdge.
     *  \details If the edge is not associated with an edge on this graph, the
     *  returned edge is null.
     *  \param e the edge to get the associated edge of.
     *  \return the associated edge. */
    CMGEdge GetEdge(const MGEdge e) const;
    
    /*! \brief Get the vertex associated with an MGVertex.
     *  \details If the vertex is not associated with a vertex of this graph,
     *  the returned vertex is null.
     *  \param v the MGVertex to get the assocaited vertex of.
     *  \return the associated vertex. */
    CMGVertex GetVertex(const MGVertex v) const;
    
    /*! \brief Get iterators across the edges of the graph.
     *  \return a pair of iterators marking the begining and end of the
     *  edges in the graph. */
    inline std::pair<EdgeIter, EdgeIter> GetEdges() const {
      return std::make_pair(_e.begin(), _e.end());
    }
    
    /*! \brief Get iterator access to the neighbours of a vertex.
     *  \details If the vertex is not a part of the graph, the range will
     *  be empty.
     *  \param v the vertex to get the neighbours of.
     *  \return a pair of iterators marking the beginning and end of the
     *  neighbours of v. */
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(const CMGVertex v) const {
      if (!HasVertex(v)) return std::make_pair(_v.end(), _v.end());
      else return std::make_pair(_n.at(v).begin(), _n.at(v).end());
    }
    
    /*! \brief Get the source vertex of an edge.
     *  \details If the edge is not a part of the graph, the returned vertex
     *  is null.
     *  \param e the edge to get the source of.
     *  \return the source vertex of the edge. */
    CMGVertex GetSource(const CMGEdge e) const;
    
    /*! \brief Get the target vertex of an edge.
     *  \details If the edge is not a part of the graph, the returned vertex
     *  is null.
     *  \param e the edge to get the target of.
     *  \return the target vertex of the edge. */
    CMGVertex GetTarget(const CMGEdge e) const;
    
    /*! \brief Get the two vertices of an edge.
     *  \details If the edge is not part of the graph, the returned vertex
     *  is null.
     *  \param e the edge to get the vertices of.
     *  \return a pair of the two vertices which the edge is between. */
    std::pair<CMGVertex, CMGVertex> GetVertices(const CMGEdge e) const;
    
    /*! \brief Get iterators across the vertices of the graph.
     *  \return a pair of iterators marking the begining and end of the
     *  vertices in the graph. */
    inline std::pair<VertIter, VertIter> GetVertices() const {
      return std::make_pair(_v.begin(), _v.end());
    }
    
    /*! \brief Check if the graph has a vertex directly associated with an
     *  MGVertex.
     *  \details Directly associated means that the vertex cannot be condensed
     *  into another vertex.
     *  \param v the vertex to check for.
     *  \return if the vertex is directly associated with the graph or not. */
    inline bool HasVertex(const MGVertex& v) const {
      if (!v) return false;
      return _vmap.find(v) != _vmap.end();
    }
    
    /*! \brief Check of the graph has a vertex associated with an MGVertex.
     *  \details Association in this case includes a vertex where one of the
     *  condensed vertices is the provided vertex.
     *  \param v the vertex to check for.
     *  \return if the vertex is associated with the graph or not. */
    bool HasCondensedVertex(const MGVertex& v) const;
    
    /*! \brief Check if the graph has a vertex.
     *  \param v the vertex to check for.
     *  \return if the vertex is part of the graph or not. */
    inline bool HasVertex(const CMGVertex v) const {
      if (!v) return false;
      return _g.HasVertex(v.get());
    }
    
    /*! \brief Check if the graph has an edge associated with an MGEdge.
     *  \param e the edge to check for.
     *  \return if the edge is associated with the graph or not. */
    inline bool HasEdge(const MGEdge& e) const {
      if (!e) return false;
      return _emap.find(e) != _emap.end();
    }
    
    /*! \brief Check if the graph has a edge.
     *  \param e the edge to check for.
     *  \return if the edge is part of the graph or not. */
    inline bool HasEdge(const CMGEdge e) const {
      if (!e) return false;
      return _g.HasEdge(e.get());
    }
    
    /*! \brief Check if the graph has an edge between two vertices.
     *  \param u, v the vertices to check for an edge between.
     *  \return if there is an edge between the two vertices or not. */
    inline bool HasEdge(const CMGVertex u, const CMGVertex v) const {
      if (!u || !v) return false;
      return _g.HasEdge(u.get(), v.get());
    }
    
    /*! \brief The number of edges in the graph.
     *  \return the number of edges. */
    inline size_ NumEdges() const { return _g.NumEdges(); }
    
    /*! \brief The number of vertices in the graph.
     *  \return the number of vertices. */
    inline size_ NumVertices() const { return _g.NumVertices(); }
    
  private:
    //! \brief Source molecule of the molecular graph.
    MolecularGraph _source;
    //! \brief Underlying graph
    graph_type _g;
    //! \brief Map MGVertex to their corresponding CMGVertex
    VertMap _vmap;
    //! \brief Map MGEdge to their corresponding CMGEdge
    EdgeMap _emap;
    //! \brief Container for giving iterator access to all vertices in graph.
    VertContain _v;
    //! \brief Container for giving iterator access to all edges in graph.
    EdgeContain _e;
    //! \brief Container for neighbours of a vertex
    NbrsContain _n;
  };
  
}

#endif  /* INDIGOX_GRAPH_CONDENSED_HPP */