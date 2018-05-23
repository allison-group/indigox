/*! \file molecular.hpp */
#ifndef INDIGOX_GRAPH_MOLECULAR_HPP
#define INDIGOX_GRAPH_MOLECULAR_HPP

#include <map>
#include <memory>
#include <set>
#include <stdexcept>

#include "base_graph.hpp"
#include "../utils/numerics.hpp"

namespace indigox{
  class IXAtom;
  class IXBond;
  class IXMolecule;
  namespace test { class IXMolecularGraph; }
  
  typedef std::shared_ptr<IXAtom> Atom;
  typedef std::shared_ptr<IXBond> Bond;
  typedef std::shared_ptr<IXMolecule> Molecule;
  
  typedef std::weak_ptr<IXAtom> _Atom;
  typedef std::weak_ptr<IXBond> _Bond;
  typedef std::weak_ptr<IXMolecule> _Molecule;
  
  namespace graph {
    
    class IXMolecularGraph;
    class IXMGVertex;
    class IXMGEdge;
    //! \brief shared_ptr for normal use of the IXMolecularGraph class.
    typedef std::shared_ptr<IXMolecularGraph> MolecularGraph;
    //! \brief shared_ptr for normal use of the IXMGVertex class.
    typedef std::shared_ptr<IXMGVertex> MGVertex;
    //! \brief shared_ptr for normal use of the IXMGEdge class.
    typedef std::shared_ptr<IXMGEdge> MGEdge;
    /*! \brief weak_ptr for non-ownership reference to the IXMolecularGraph
     *  class.
     *  \details Intended for internal use only. */
    typedef std::weak_ptr<IXMolecularGraph> _MolecularGraph;
    /*! \brief weak_ptr for non-ownership reference to the IXMGVertex class.
     *  \details Intended for internal use only. */
    typedef std::weak_ptr<IXMGVertex> _MGVertex;
    /*! \brief weak_ptr for non-ownership reference to the IXMGEdge  .
     *  \details Intended for internal use only. */
    typedef std::weak_ptr<IXMGEdge> _MGEdge;
    
    /*! \brief Class for the vertices of a IXMolecularGraph. */
    class IXMGVertex : public std::enable_shared_from_this<IXMGVertex> {
      //! \brief Friendship allows IXMolecularGraph to add vertices.
      friend class IXMolecularGraph;
    public:
      IXMGVertex() = delete;  // no default constructor
      /*! \brief Get the atom associated with this vertex.
       *  \return the atom associated with this vertex, if it is still alive. */
      Atom GetAtom() const { return _atom.lock(); }
    private:
      /*! \brief Construct a vertex from an atom.
       *  \details Private constructor ensures that only IXMolecularGraph can
       *  create IXMGVertex instances.
       *  \param a the atom to associate with this vertex. */
      IXMGVertex(Atom a) : _atom(a) { }
      
      //! \brief Reference to the atom this edge is associated with.
      _Atom _atom;
    };
    
    /*! \brief Class for the edges of a IXMolecularGraph. */
    class IXMGEdge : public std::enable_shared_from_this<IXMGEdge> {
      //! \brief Friendship allows IXMolecularGraph to add edges.
      friend class IXMolecularGraph;
    public:
      IXMGEdge() = delete;  // no default constructor
      /*! \brief Get the bond associated with this edge.
       *  \return the bond associated with this edge, if it is still alive. */
      Bond GetBond() const { return _bond.lock(); }
    private:
      /*! \brief Construct an edge from a bond.
       *  \details Private constructor ensures that only IXMolecularGraph can
       *  create IXMGEdge instances.
       *  \param b the bond to associate with this edge. */
      IXMGEdge(Bond b) : _bond(b) { }
      
      //! \brief Reference to the bond this edge is associated with.
      _Bond _bond;
    };
    
    /*! \brief Class containing a graph representation of a molecule.
     *  \details The IXMolecularGraph is designed to be maintained by the
     *  IXMolecule instance owning it. To that end, all the modifying methods
     *  assume that the parameters feed to them are valid. However, all the
     *  accessing methods do not make this assumption and so perform sanity
     *  checks. */
    class IXMolecularGraph {
      //! \brief Friendship allows an IXMolecule to own a graph.
      friend class indigox::IXMolecule;
      //! \brief Friendship allows IXMolecularGraph to be tested.
      friend class indigox::test::IXMolecularGraph;
      //! \brief Type of the underlying IXGraphBase
      typedef IXGraphBase<IXMGVertex, IXMGEdge> graph_type;
      //! \brief Type of the iterator returned by GetEdges() method.
      typedef std::vector<MGEdge>::const_iterator EdgeIter;
      //! \brief Type of the iterator returned by GetVertices() method.
      typedef std::vector<MGVertex>::const_iterator VertIter;
      //! \brief Type of the iterator returned by GetNeighbours() method.
      typedef std::vector<MGVertex>::const_iterator NbrsIter;
      //! \brief Type of the iterator over components of the graph
      typedef std::vector<std::vector<MGVertex>>::const_iterator CompIter;
      
    public:
      IXMolecularGraph() = delete;  // no default constructor
      
    private:
      /*! \brief Construct with a molecule.
       *  \param mol the molecule to reference to. */
      IXMolecularGraph(const Molecule mol) : _source(mol), _g() { }
      
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
       *  \details Iterators are valid for the graph state at the point the
       *  method was called. Subsequent method calls will invalidate all
       *  previous iterators.
       *  \return a pair of iterators marking the begining and end of the
       *  edges in the graph. */
      std::pair<EdgeIter, EdgeIter> GetEdges();
      
      /*! \brief Get iterator access to the neighbours of a vertex.
       *  \details Iterators are valid for the graph state at the point the
       *  method was called. Any subsequent method calls will invalidate all
       *  previously returned iterators, regardless of what vertex they were
       *  called for. If the vertex is not a part of the graph, the range will
       *  be empty.
       *  \param v the vertex to get the neighbours of.
       *  \return a pair of iterators marking the beginning and end of the
       *  neighbours of v. */
      std::pair<NbrsIter, NbrsIter> GetNeighbours(const MGVertex v);
      
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
       *  \details Iterators are valid for the graph state at the point the
       *  method was called. Subsequent method calls  will invalidate all
       *  previous iterators.
       *  \return a pair of iterators marking the begining and end of the
       *  vertices in the graph. */
      std::pair<VertIter, VertIter> GetVertices();
      
      /*! \brief Check if the graph has a vertex associated with an atom.
       *  \param v that atom to check for.
       *  \return if the atom is associated with the graph or not. */
      inline bool HasVertex(const Atom v) const {
        if (!v) return false;
        return _verts.find(v) != _verts.end();
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
        return _edges.find(e) != _edges.end();
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
      
      /*! \brief Determine if the graph is connected.
       *  \details A graph is connected if it contains only one connected
       *  component. Here, a graph is also deemed to be connected if it contains
       *  no vertices.
       *  \return if the graph is connected or not. */
      inline bool IsConnected() { return _g.NumConnectedComponents() <= 1; }
      
      /*! \brief Calculate number of connected components.
       *  \return the number of connected components of the graph. */
      inline size_ NumConnectedComponents() { return _g.NumConnectedComponents(); }
      
      /*! \brief Calculate the connected components.
       *  \return a pair of iterators providing access to the components. Each
       *  component is a std::vector<MGVertex>. */
      std::pair<CompIter, CompIter> GetConnectedComponents();
      
    private:
      //! \brief Source molecule of the molecular graph.
      _Molecule _source;
      //! \brief Underlying graph
      graph_type _g;
      //! \brief Map Atoms to their corresponding MGVertex
      std::map<Atom, MGVertex> _verts;
      //! \brief Map Bonds to their corresponding MGEdge
      std::map<Bond, MGEdge> _edges;
      //! \brief Container for giving iterator access to all vertices in graph.
      std::vector<MGVertex> _vert_access;
      //! \brief Container for giving iterator access to all edges in graph.
      std::vector<MGEdge> _edge_access;
      //! \brief Container for neighbours of a vertex
      std::vector<MGVertex> _nbrs_access;
      //! \brief Container for components of the graph
      std::vector<std::vector<MGVertex>> _components;
    
    };
  }  // namespace graph
}  // namespace indigox

#endif /* INDIGOX_GRAPH_MOLECULAR_HPP */
