/*! \file molecular.hpp */
#ifndef INDIGOX_GRAPH_MOLECULAR_HPP
#define INDIGOX_GRAPH_MOLECULAR_HPP

#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

#include <EASTL/vector_map.h>

#include "base_graph.hpp"
#include "../algorithm/graph/cycles.hpp"
#include "../utils/common.hpp"
#include "../utils/numerics.hpp"

// Forward declares
namespace indigox {
  class IXAtom;
  class IXBond;
  class IXMolecule;
  namespace test {
    struct TestMolecularGraph;
    struct TestMolecularVertex;
    struct TestMolecularEdge;
  }
  
  using Atom = std::shared_ptr<IXAtom>;
  using Bond = std::shared_ptr<IXBond>;
  using Molecule = std::shared_ptr<IXMolecule>;
  
  using _Atom = std::weak_ptr<IXAtom>;
  using _Bond = std::weak_ptr<IXBond>;
  using _Molecule = std::weak_ptr<IXMolecule>;
  
  namespace graph {
    class IXMolecularGraph;
    class IXMGVertex;
    class IXMGEdge;
    
    //! \brief shared_ptr for normal use of the IXMolecularGraph class.
    using MolecularGraph = std::shared_ptr<IXMolecularGraph>;
    //! \brief shared_ptr for normal use of the IXMGVertex class.
    using MGVertex = std::shared_ptr<IXMGVertex>;
    //! \brief shared_ptr for normal use of the IXMGEdge class.
    using MGEdge = std::shared_ptr<IXMGEdge>;
    /*! \brief weak_ptr for non-ownership reference to the IXMolecularGraph
     *  class.
     *  \details Intended for internal use only. */
    using _MolecularGraph = std::weak_ptr<IXMolecularGraph>;
    /*! \brief weak_ptr for non-ownership reference to the IXMGVertex class.
     *  \details Intended for internal use only. */
    using _MGVertex = std::weak_ptr<IXMGVertex>;
    /*! \brief weak_ptr for non-ownership reference to the IXMGEdge  .
     *  \details Intended for internal use only. */
    using _MGEdge = std::weak_ptr<IXMGEdge>;
  }
}

// Local declarations
namespace indigox::graph {
  /*! \brief Class for the vertices of a IXMolecularGraph. */
  class IXMGVertex : public std::enable_shared_from_this<IXMGVertex> {
    //! \brief Friendship allows IXMolecularGraph to add vertices.
    friend class IXMolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;
    //! \brief Friendship allows for testing
    friend struct indigox::test::TestMolecularVertex;
    
  public:
    /*! \brief Get the atom associated with this vertex.
     *  \return the atom associated with this vertex, if it is still alive. */
    inline Atom GetAtom() const { return _atom.lock(); }
    
    /*! \brief Get the graph this vertex is part of.
     *  \return the owning graph. */
    inline MolecularGraph GetGraph() const { return _graph.lock(); }
    
    IXMGVertex() = delete;  // no default constructor
    
  private:
    /*! \brief Construct a vertex from an atom.
     *  \details Private constructor ensures that only IXMolecularGraph can
     *  create IXMGVertex instances.
     *  \param a the atom to associate with this vertex. */
    IXMGVertex(Atom a, MolecularGraph graph) : _atom(a), _graph(graph) { }
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXMGVertex>& construct,
                                   const uint32_t version);
    
    //! \brief Reference to the atom this edge is associated with.
    _Atom _atom;
    //! \brief Reference to the molecular graph
    _MolecularGraph _graph;
  };
  
  /*! \brief Class for the edges of a IXMolecularGraph. */
  class IXMGEdge : public std::enable_shared_from_this<IXMGEdge> {
    //! \brief Friendship allows IXMolecularGraph to add edges.
    friend class IXMolecularGraph;
    //! \brief Friendship allows for serialisation
    friend class cereal::access;
    //! \brief Friendship allows for testing
    friend struct indigox::test::TestMolecularEdge;
    
  public:
    /*! \brief Get the bond associated with this edge.
     *  \return the bond associated with this edge, if it is still alive. */
    Bond GetBond() const { return _bond.lock(); }
    
    /*! \brief Get the graph this edge is part of.
     *  \return the owning graph. */
    inline MolecularGraph GetGraph() const { return _graph.lock(); }
    
    IXMGEdge() = delete;  // no default constructor
    
  private:
    /*! \brief Construct an edge from a bond.
     *  \details Private constructor ensures that only IXMolecularGraph can
     *  create IXMGEdge instances.
     *  \param b the bond to associate with this edge. */
    IXMGEdge(Bond b, MolecularGraph graph) : _bond(b), _graph(graph) { }
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXMGEdge>& construct,
                                   const uint32_t version);
    
    //! \brief Reference to the bond this edge is associated with.
    _Bond _bond;
    //! \brief Reference to the molecular graph
    _MolecularGraph _graph;
  };
  
  /*! \brief Class containing a graph representation of a molecule.
   *  \details The IXMolecularGraph is designed to be maintained by the
   *  IXMolecule instance owning it. To that end, all the modifying methods
   *  assume that the parameters feed to them are valid. However, all the
   *  accessing methods do not make this assumption and so perform sanity
   *  checks. */
  class IXMolecularGraph
  : public std::enable_shared_from_this<IXMolecularGraph> {
  public:
    //! \brief Friendship allows an IXMolecule to own a graph.
    friend class indigox::IXMolecule;
    //! \brief Friendship allows IXMolecularGraph to be tested.
    friend struct indigox::test::TestMolecularGraph;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    //! \brief Friendship allows graph algorithms access to underlying graph
    friend struct indigox::graph::access;
    
    //! \brief Type of the underlying IXGraphBase
    using graph_type = IXGraphBase<IXMGVertex, IXMGEdge>;
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
    
  private:
    /*! \brief Construct with a molecule.
     *  \param mol the molecule to reference to. */
    IXMolecularGraph(const Molecule mol) : _source(mol), _g() { }
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXMolecularGraph>& construct,
                                   const uint32_t version);
    
    // modifcation methods are private so that the structure of the graph can
    // be controlled only by the molecule owning it.
    /*! \brief Add an edge to the graph.
     *  \details If either of the atoms of the bond are not part of the graph
     *  already, they are also added to it. It is assumed that the provided
     *  bond is not already associated with another edge.
     *  \param bnd the bond the new edge is associated with.
     *  \return shared_ptr to the newly added edge. */
    MGEdge AddEdge(const Bond bnd);
    
    /*! \brief Add a vertex to the graph.
     *  \details It is assumed that the provided atom is not already
     *  associated with another vertex.
     *  \param atm the atom the new vertex is associated with.
     *  \return shared_ptr to the newly added vertex. */
    MGVertex AddVertex(const Atom atm);
    
    //! \brief Clears all vertices and edges from the graph.
    void Clear();
    
    /*! \brief Remove an edge from the graph.
     *  \details It is assumed that the provided edge is a part of the graph.
     *  \param e the edge to remove. */
    void RemoveEdge(const MGEdge e);
    
    /*! \brief Remove an edge from between two vertices.
     *  \details It is assumed that there is an edge between the provided
     *  vertices.
     *  \param u, v the vertices to remove the edge from between. */
    void RemoveEdge(const MGVertex u, const MGVertex v);
    
    /*! \brief Remove a vertex from the graph.
     *  \details Removing a vertex also removes all edges it is involved in.
     *  It is assumed that the provided vertex is part of the graph.
     *  \param v the vertex to remove. */
    void RemoveVertex(const MGVertex v);
    
  public:
    IXMolecularGraph() = default;  // no default constructor
    
    /*! \brief Create a snapshot of the graph in its current state.
     *  \details No link to a molecule is maintained.
     *  \return snapshot of the graph. */
    MolecularGraph CreateSnapshot() const;
    
    /*! \brief Induce a subgraph from the range of vertices.
     *  \details Induced subgraph has the same vertices and edges as its parent
     *  graph. Additionally, its source Molecule is the same. This is a
     *  vertex induced subgraph, meaning that all edges where both vertices are
     *  in the provided range will be in the induced graph.
     *  \tparam InputIt type of the iterator range provided.
     *  \param begin,end marking the range of vertices to induce subgraph on.
     *  \return a new MolecularGraph. */
    template <class InputIt>
    MolecularGraph InduceSubgraph(InputIt begin, InputIt end) const {
      MolecularGraph G = std::make_shared<IXMolecularGraph>();
      G->_source = _source;
      for (auto& vs : _at2v) {
        if (std::find(begin, end, vs.second) == end) continue;
        G->_g.AddVertex(vs.second.get());
        G->_at2v.emplace(vs.first, vs.second);
        G->_v.emplace_back(vs.second);
        G->_n.emplace(vs.second, NbrsContain::mapped_type());
      }
      
      for (auto& es : _bn2e) {
        MGVertex u = GetSource(es.second);
        MGVertex v = GetTarget(es.second);
        if (!G->HasVertex(u) || !G->HasVertex(v)) continue;
        G->_g.AddEdge(u.get(), v.get(), es.second.get());
        G->_bn2e.emplace(es.first, es.second);
        G->_e.emplace_back(es.second);
        G->_n[u].emplace_back(v);
        G->_n[v].emplace_back(u);
      }
      
      return G;
    }
    
    /*! \brief The degree of a vertex.
     *  \details If the vertex is not part of the graph, the returned value is
     *  std::numeric_limits<size_>::max().
     *  \param v the vertex to obtain the degree of.
     *  \return the degree of the vertex. */
    size_ Degree(const MGVertex v) const;
    
    /*! \brief Get the edge between two atoms.
     *  \details If there is no edge between the vertices, or at least one of
     *  the vertices is not part of the graph, the returned edge is null.
     *  \param u, v the vertices to get the edge between.
     *  \return the edge between the two vertices. */
    MGEdge GetEdge(const MGVertex u, const MGVertex v) const;
    
    /*! \brief Get the edge associated with a bond.
     *  \details If the bond is not associated with an edge on this graph, the
     *  returned edge is null.
     *  \param bnd the bond to get the associated edge of.
     *  \return the associated edge. */
    MGEdge GetEdge(const Bond bnd) const;
    
    /*! \brief Get the vertex associated with an atom.
     *  \details If the atom is not associated with a vertex of this graph,
     *  the returned vertex is null.
     *  \param atm the atom to get the assocaited vertex of.
     *  \return the associated vertex. */
    MGVertex GetVertex(const Atom atm) const;
    
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
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(const MGVertex v) const {
      if (!HasVertex(v)) return std::make_pair(_v.end(), _v.end());
      else return std::make_pair(_n.at(v).begin(), _n.at(v).end());
    }
    
    /*! \brief Get the source vertex of an edge.
     *  \details If the edge is not a part of the graph, the returned vertex
     *  is null.
     *  \param e the edge to get the source of.
     *  \return the source vertex of the edge. */
    MGVertex GetSource(const MGEdge e) const;
    
    /*! \brief Get the target vertex of an edge.
     *  \details If the edge is not a part of the graph, the returned vertex
     *  is null.
     *  \param e the edge to get the target of.
     *  \return the target vertex of the edge. */
    MGVertex GetTarget(const MGEdge e) const;
    
    /*! \brief Get the two vertices of an edge.
     *  \details If the edge is not part of the graph, the returned vertex
     *  is null.
     *  \param e the edge to get the vertices of.
     *  \return a pair of the two vertices which the edge is between. */
    std::pair<MGVertex, MGVertex> GetVertices(const MGEdge e) const;
    
    /*! \brief Get iterators across the vertices of the graph.
     *  \return a pair of iterators marking the begining and end of the
     *  vertices in the graph. */
    inline std::pair<VertIter, VertIter> GetVertices() const {
      return std::make_pair(_v.begin(), _v.end());
    }
    
    /*! \brief Check if the graph has a vertex associated with an atom.
     *  \param v that atom to check for.
     *  \return if the atom is associated with the graph or not. */
    inline bool HasVertex(const Atom v) const {
      if (!v) return false;
      return _at2v.find(v) != _at2v.end();
    }
    
    /*! \brief Check if the graph has a vertex.
     *  \param v the vertex to check for.
     *  \return if the vertex is part of the graph or not. */
    inline bool HasVertex(const MGVertex v) const {
      if (!v) return false;
      return _g.HasVertex(v.get());
    }
    
    /*! \brief Check if the graph has an edge associated with a bond.
     *  \param e the bond to check for.
     *  \return if the bond is associated with the graph or not. */
    inline bool HasEdge(const Bond e) const {
      if (!e) return false;
      return _bn2e.find(e) != _bn2e.end();
    }
    
    /*! \brief Check if the graph has a edge.
     *  \param e the edge to check for.
     *  \return if the edge is part of the graph or not. */
    inline bool HasEdge(const MGEdge e) const {
      if (!e) return false;
      return _g.HasEdge(e.get());
    }
    
    /*! \brief Check if the graph has an edge between two vertices.
     *  \param u, v the vertices to check for an edge between.
     *  \return if there is an edge between the two vertices or not. */
    inline bool HasEdge(const MGVertex u, const MGVertex v) const {
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
    _Molecule _source;
    //! \brief Underlying graph
    graph_type _g;
    //! \brief Map Atoms to their corresponding MGVertex
    AtomMap _at2v;
    //! \brief Map Bonds to their corresponding MGEdge
    BondMap _bn2e;
    //! \brief Container for giving iterator access to all vertices in graph.
    VertContain _v;
    //! \brief Container for giving iterator access to all edges in graph.
    EdgeContain _e;
    //! \brief Container for neighbours of a vertex
    NbrsContain _n;
  };
}  // namespace indigox::graph

#endif /* INDIGOX_GRAPH_MOLECULAR_HPP */
