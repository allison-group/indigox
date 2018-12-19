/*! \file condensed.hpp */
#ifndef INDIGOX_GRAPH_CONDENSED_HPP
#define INDIGOX_GRAPH_CONDENSED_HPP

#include "../utils/fwd_declares.hpp"
#include "base_graph.hpp"

#include <EASTL/bitset.h>
#include <EASTL/vector_map.h>
#include <EASTL/vector_set.h>
#include <memory>
#include <vector>

namespace indigox::graph {
  //! \brief type used to store the isomorphism testing mask for IXCMGVertex.
  using VertexIsoMask = eastl::bitset<64, uint64_t>;
  //! \brief type used to store the isomorphism testing mask for IXCMGEdge.
  using EdgeIsoMask = eastl::bitset<16, uint16_t>;

  /*! \brief Class for the vertices of an IXCondensedMolecularGraph. */
  class CMGVertex {
  public:
    //! \brief Enum for the contracted symmetry groups
    enum class ContractedSymmetry {
      Hydrogen, //!< Group for contracted hydrogen atoms.
      Fluorine, //!< Group for contracted fluorine atoms.
      Chlorine, //!< Group for contracted chlorine atoms.
      Bromine,  //!< Group for contracted bromine atoms.
      Iodine    //!< Group for contracted iodine atoms.
    };

  private:
    //! \brief Friendship allows IXCondensedMolecularGraph to add vertices.
    friend class CondensedMolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;
    //! \brief Type used to store condensed vertices
    using CondensedVertex = std::pair<ContractedSymmetry, MGVertex>;

    struct CMGVertexData;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(CMGVertex);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(CMGVertex, v);

  private:
    /*! \brief Construct a vertex from an MGVertex.
     *  \param v the MGVertex to associate with this vertex.
     *  \param g the condensed graph this vertex will be a part of. */
    CMGVertex(const MGVertex &v, const CondensedMolecularGraph &g);

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);


  public:
    /*! \brief Get the MGVertex associated with this vertex.
     *  \return the associated MGVertex, if it is still alive. */
    const MGVertex& GetSource() const;

    /*! \brief Get the graph this vertex is part of.
     *  \return the owning graph. */
    const CondensedMolecularGraph &GetGraph() const;

    /*! \brief Get the number of contracted vertices.
     *  \return the number of contracted vertices. */
    size_t NumContracted() const;

    /*! \brief Get the number of contracted vertices in the symmetry group.
     *  \return the number of vertices of the given symmetry group contracted.
     */
    size_t NumContracted(ContractedSymmetry sym) const;

    /*! \brief Get the isomorphism testing mask.
     *  \return the isomorphism testing mask. */
    const VertexIsoMask &GetIsomorphismMask() const;

    /*! \brief Checks if a given MGVertex is contracted into this vertex.
     *  \param v the vertex to check for.
     *   return if the provided vertex is contracted into this one. */
    bool IsContractedHere(const MGVertex &v) const;

    std::vector<MGVertex> GetContractedVertices() const;
    const std::vector<CondensedVertex> &GetCondensedVertices() const;

  private:
    //! \brief The vertex data.
    std::shared_ptr<CMGVertexData> m_data;
  };
  using ContractedSymmetry = CMGVertex::ContractedSymmetry;

  /*! \brief Class for the edges of an IXCondensedMolecularGraph. */
  class CMGEdge {
    //! \brief Friendship allows IXCondensedMolecularGraph to add edges.
    friend class CondensedMolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;

    struct CMGEdgeData;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(CMGEdge);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(CMGEdge, e);

  private:
    /*! \brief Construct an edge from an MGEdge.
     *  \param e the MGEdge to associate with this edge.
     *  \param g the condensed graph this edge will be a part of. */
    CMGEdge(const MGEdge &e, const CondensedMolecularGraph &g);

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

  public:
    /*! \brief Get the MGEdge associated with this vertex.
     *  \return the associated MGEdge, if it is atill alive. */
    const MGEdge& GetSource() const;

    /*! \brief Get the graph this edge is part of.
     *  \return the owning graph. */
    const CondensedMolecularGraph &GetGraph() const;

    /*! \brief Get the isomorphism testing mask.
     *  \return the isomorphism testing mask. */
    const EdgeIsoMask &GetIsomorphismMask() const;

  private:
    //! \brief Source edge
    std::shared_ptr<CMGEdgeData> m_data;
  };

  class CondensedMolecularGraph
      : public BaseGraph<CMGVertex, CMGEdge, CondensedMolecularGraph> {
  public:
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    //! \brief Friendship allows for generating from a source
    friend CondensedMolecularGraph Condense(const MolecularGraph &);
    friend class MolecularGraph;

    //! \brief Type of the underlying IXGraphBase
    using graph_type = BaseGraph<CMGVertex, CMGEdge, CondensedMolecularGraph>;
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
    //! \brief Type used for vertices
    using VertexType = CMGVertex;
    //! \brief Type used for edges
    using EdgeType = CMGEdge;

  private:
    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

    // Modification methods are private so that the CMG is a snapshot of the MG
    // at time of creation
    /*! \brief Add an edge to the graph.
     *  \details Adds an edge between previously added vertices. Assumes that
     *  the source and target vertex of the provided MGEdge have already been
     *  added to the graph.
     *  \param e the source MGEdge.
     *  \return shared_ptr tp the newly added edge. */
    CMGEdge AddEdge(const MGEdge &e);

    /*! \brief Add a vertex to the graph.
     *  \details It is assumed that the provided source vertex is not viable
     *  for condensing.
     *  \param v the source MGVertex.
     *  \return shared_ptr to the newly added vertex. */
    CMGVertex AddVertex(const MGVertex &v);

//    void Clear();

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(CondensedMolecularGraph);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(CondensedMolecularGraph, G);

  private:
    /*! \brief Construct a condensed molecular graph from a MolecularGraph.
     *  \param g the molecular graph to construct from. */
    CondensedMolecularGraph(const MolecularGraph &g);

  public:
    /*! \brief Induce a subgraph from the range of vertices.
     *  \details Induced subgraph has the same vertices and edges as its parent
     *  graph. Additionally, its source MolecularGraph is the same. This is a
     *  vertex induced subgraph, meaning that all edges where both vertices are
     *  in the provided range will be in the induced graph.
     *  \tparam InputIt type of the iterator range provided.
     *  \param begin,end marking the range of vertices to induce subgraph on.
     *  \return a new CondensedMolecularGraph. */
    CondensedMolecularGraph Subgraph(std::vector<CMGVertex> &verts);

    /*! \brief Create a subgraph from the range of vertices and edges.
     *  \details Subgraph has the same vertices and edges as its parent graph.
     *  Additionally, its source MolecularGraph is the same. This subgraph is
     *  such that only edges within the provided range are added to it, as long
     *  as both vertices are in the provided vertex range.
     *  \tparam VertIt type of the vertex iterator range provided.
     *  \tparam EdgeIt type of the edge iterator range provided.
     *  \param v_begin,v_end marking the range of vertices to create subgraph.
     *  \param e_begin,e_end marking the range of edges to include in subgraph.
     *  \return a new CondensedMolecularGraph. */
    CondensedMolecularGraph Subgraph(std::vector<CMGVertex> &verts,
                                      std::vector<CMGEdge> &edges);

    bool IsSubgraph() const;
//     {
//      return bool(_subg);
//    }

    /*! \brief Get the source MolecularGraph.
     *  \return the molecular graph used to construt this. */
    const MolecularGraph &GetMolecularGraph() const;

    const CondensedMolecularGraph &GetSuperGraph() const;

    using graph_type::GetEdge;
    using graph_type::HasEdge;
    using graph_type::HasVertex;

    /*! \brief Get the edge associated with an MGEdge.
     *  \details If the edge is not associated with an edge on this graph, the
     *  returned edge is null.
     *  \param e the edge to get the associated edge of.
     *  \return the associated edge. */
    const CMGEdge& GetEdge(const MGEdge &e) const;

    /*! \brief Get the vertex associated with an MGVertex.
     *  \details If the vertex is not associated with a vertex of this graph,
     *  the returned vertex is null.
     *  \param v the MGVertex to get the assocaited vertex of.
     *  \return the associated vertex. */
    const CMGVertex& GetVertex(const MGVertex &v) const;
    const CMGVertex& GetCondensedVertex(const MGVertex &v) const;

    /*! \brief Check if the graph has a vertex directly associated with an
     *  MGVertex.
     *  \details Directly associated means that the vertex cannot be condensed
     *  into another vertex.
     *  \param v the vertex to check for.
     *  \return if the vertex is directly associated with the graph or not. */
    bool HasVertex(const MGVertex &v) const;

    /*! \brief Check of the graph has a vertex associated with an MGVertex.
     *  \details Association in this case includes a vertex where one of the
     *  condensed vertices is the provided vertex.
     *  \param v the vertex to check for.
     *  \return if the vertex is associated with the graph or not. */
    bool HasCondensedVertex(const MGVertex &v) const;

    /*! \brief Check if the graph has an edge associated with an MGEdge.
     *  \param e the edge to check for.
     *  \return if the edge is associated with the graph or not. */
    bool HasEdge(const MGEdge &e) const;

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
  };

  CondensedMolecularGraph Condense(const MolecularGraph &G);

} // namespace indigox::graph

#endif /* INDIGOX_GRAPH_CONDENSED_HPP */
