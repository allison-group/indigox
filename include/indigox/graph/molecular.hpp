/*! \file molecular.hpp */
#ifndef INDIGOX_GRAPH_MOLECULAR_HPP
#define INDIGOX_GRAPH_MOLECULAR_HPP

#include "../algorithm/graph/cycles.hpp"
#include "../utils/fwd_declares.hpp"
#include "../utils/modifable_object.hpp"
#include "base_graph.hpp"

#include <EASTL/vector_map.h>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

// Local declarations
namespace indigox::graph {
  /*! \brief Class for the vertices of a IXMolecularGraph. */
  class MGVertex {
    //! \brief Friendship allows IXMolecularGraph to add vertices.
    friend class MolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;

  public:
    /*! \brief Get the atom associated with this vertex.
     *  \return the atom associated with this vertex, if it is still alive. */
    const Atom &GetAtom() const;

    /*! \brief Get the graph this vertex is part of.
     *  \return the owning graph. */
    const MolecularGraph &GetGraph() const;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(MGVertex);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(MGVertex, v);

  private:
    /*! \brief Construct a vertex from an atom.
     *  \details Private constructor ensures that only MolecularGraph can
     *  create valid MGVertex instances.
     *  \param a the atom to associate with this vertex. */
    MGVertex(const Atom &a, const MolecularGraph &graph);

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

  private:
    struct MGVertexData;
    std::shared_ptr<MGVertexData> m_data;
  };

  /*! \brief Class for the edges of a IXMolecularGraph. */
  class MGEdge {
    //! \brief Friendship allows IXMolecularGraph to add edges.
    friend class MolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;

  public:
    /*! \brief Get the bond associated with this edge.
     *  \return the bond associated with this edge, if it is still alive. */
    const Bond& GetBond() const;

    /*! \brief Get the graph this edge is part of.
     *  \return the owning graph. */
    const MolecularGraph &GetGraph() const;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(MGEdge);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(MGEdge, e);

  private:
    /*! \brief Construct an edge from a bond.
     *  \details Private constructor ensures that only IXMolecularGraph can
     *  create IXMGEdge instances.
     *  \param b the bond to associate with this edge. */
    MGEdge(const Bond &b, const MolecularGraph &graph);

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

  private:
    struct MGEdgeData;
    std::shared_ptr<MGEdgeData> m_data;
  };

  /*! \brief Class containing a graph representation of a molecule.
   *  \details The IXMolecularGraph is designed to be maintained by the
   *  IXMolecule instance owning it. To that end, all the modifying methods
   *  assume that the parameters feed to them are valid. However, all the
   *  accessing methods do not make this assumption and so perform sanity
   *  checks. */
  class MolecularGraph : public BaseGraph<MGVertex, MGEdge, MolecularGraph>{
  public:
    //! \brief Friendship allows an Molecule to own a graph.
    friend class indigox::Molecule;
    //! \brief Friendship allows serialisation
    friend class cereal::access;

    //! \brief Type of the underlying IXGraphBase
    using graph_type = BaseGraph<MGVertex, MGEdge, MolecularGraph>;
    //! \brief Container for vertices
    using VertContain = std::vector<MGVertex>;
    //! \brief Container for edges
    using EdgeContain = std::vector<MGEdge>;
    //! \brief Container for neighbours of vertices
    using NbrsContain = std::map<MGVertex, VertContain>;
    //! \brief Container for mapping atoms to vertices
    using AtomMap = eastl::vector_map<Atom, MGVertex>;
    //! \brief Container for mapping bonds to edges
    using BondMap = eastl::vector_map<Bond, MGEdge>;

  public:
    //! \brief Type of the iterator returned by GetEdges() method.
    using EdgeIter = EdgeContain::const_iterator;
    //! \brief Type of the iterator returned by GetVertices() method.
    using VertIter = VertContain::const_iterator;
    //! \brief Type of the iterator returned by GetNeighbours() method.
    using NbrsIter = NbrsContain::mapped_type::const_iterator;
    //! \brief Type used for vertices
    using VertexType = MGVertex;
    //! \brief Type used for edges
    using EdgeType = MGEdge;

  public:
    /*! \brief Construct with a molecule.
     *  \param mol the molecule to reference to. */
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(MolecularGraph);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(MolecularGraph, g);

  public:
    MolecularGraph(const Molecule &mol);

  private:
    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

    // modifcation methods are private so that the structure of the graph can
    // be controlled only by the molecule owning it.
    /*! \brief Add an edge to the graph.
     *  \details If either of the atoms of the bond are not part of the graph
     *  already, they are also added to it. It is assumed that the provided
     *  bond is not already associated with another edge.
     *  \param bnd the bond the new edge is associated with.
     *  \return shared_ptr to the newly added edge. */
    MGEdge AddEdge(const Bond &bnd);

    /*! \brief Add a vertex to the graph.
     *  \details It is assumed that the provided atom is not already
     *  associated with another vertex.
     *  \param atm the atom the new vertex is associated with.
     *  \return shared_ptr to the newly added vertex. */
    MGVertex AddVertex(const Atom &atm);

    /*! \brief Remove an edge from the graph.
     *  \details It is assumed that the provided edge is a part of the graph.
     *  \param e the edge to remove. */
    void RemoveEdge(const MGEdge &e);

    /*! \brief Remove an edge from between two vertices.
     *  \details It is assumed that there is an edge between the provided
     *  vertices.
     *  \param u, v the vertices to remove the edge from between. */
    void RemoveEdge(const MGVertex &u, const MGVertex &v);

    /*! \brief Remove a vertex from the graph.
     *  \details Removing a vertex also removes all edges it is involved in.
     *  It is assumed that the provided vertex is part of the graph.
     *  \param v the vertex to remove. */
    void RemoveVertex(const MGVertex &v);

  public:
    /*! \brief Induce a subgraph from the range of vertices.
     *  \details Induced subgraph has the same vertices and edges as its parent
     *  graph. Additionally, its source Molecule is the same. This is a
     *  vertex induced subgraph, meaning that all edges where both vertices are
     *  in the provided range will be in the induced graph.
     *  \tparam InputIt type of the iterator range provided.
     *  \param begin,end marking the range of vertices to induce subgraph on.
     *  \return a new MolecularGraph. */
    MolecularGraph Subgraph(std::vector<MGVertex> &vertices);

    MolecularGraph Subgraph(std::vector<MGVertex> &vertices,
                             std::vector<MGEdge> &edges);

    bool IsSubgraph() const;
//    {
//      return !_subg.expired();
//    }

    using graph_type::GetEdge;
    using graph_type::HasEdge;
    using graph_type::HasVertex;

    /*! \brief Get the edge associated with a bond.
     *  \details If the bond is not associated with an edge on this graph, the
     *  returned edge is null.
     *  \param bnd the bond to get the associated edge of.
     *  \return the associated edge. */
    const MGEdge& GetEdge(const Bond &bnd) const;

    /*! \brief Get the vertex associated with an atom.
     *  \details If the atom is not associated with a vertex of this graph,
     *  the returned vertex is null.
     *  \param atm the atom to get the assocaited vertex of.
     *  \return the associated vertex. */
    const MGVertex& GetVertex(const Atom &atm) const;

    /*! \brief Check if the graph has a vertex associated with an atom.
     *  \param v that atom to check for.
     *  \return if the atom is associated with the graph or not. */
    bool HasVertex(const Atom &v) const;

    /*! \brief Check if the graph has an edge associated with a bond.
     *  \param e the bond to check for.
     *  \return if the bond is associated with the graph or not. */
    bool HasEdge(const Bond &e) const;

    const Molecule &GetMolecule() const;

    const MolecularGraph &GetSuperGraph() const;

    // Can only generate/get the condensed graph when molecule is frozen
    const CondensedMolecularGraph &GetCondensedGraph();

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
    
//    //! \brief Map Atoms to their corresponding MGVertex
//    AtomMap _at2v;
//    //! \brief Map Bonds to their corresponding MGEdge
//    BondMap _bn2e;
//    //! \brief Source molecule
//    Molecule _mol;
//    //! \brief Condensed version of graph
//    sCondensedMolecularGraph _cond;
//    //! \brief If is subgraph
//    wMolecularGraph _subg;
  };
} // namespace indigox::graph

#endif /* INDIGOX_GRAPH_MOLECULAR_HPP */
