/*! \file molecular.hpp */
#ifndef INDIGOX_GRAPH_MOLECULAR_HPP
#define INDIGOX_GRAPH_MOLECULAR_HPP

#include <stdexcept>

#include "base.hpp"
#include "../utils/counter.hpp"
#include "../utils/numerics.hpp"

//! Base indigoX namespace.
namespace indigox{
  class IXAtom;
  class IXBond;
  class IXMolecule;
  typedef std::shared_ptr<IXAtom> Atom;
  typedef std::shared_ptr<IXBond> Bond;
  typedef std::shared_ptr<IXMolecule> Molecule;
  typedef std::weak_ptr<IXAtom> _Atom;
  typedef std::weak_ptr<IXBond> _Bond;
  typedef std::weak_ptr<IXMolecule> _Molecule;
  
  //! Namespace for all graph related functionality.
  namespace graph {
    
    class IXMolecularGraph;
    //! \brief shared_ptr for normal use of the IXMolecularGraph class.
    typedef std::shared_ptr<IXMolecularGraph> MolecularGraph;
    /*! \brief weak_ptr for non-ownership reference to the IXMolecularGraph
     *  class.
     *  \details Intended for internal use only. */
    typedef std::weak_ptr<IXMolecularGraph> _MolecularGraph;
    
    class MGVertex {
    public:
      MGVertex() : _atom(), _id(++_id_generate) { }
      MGVertex(Atom a) : _atom(a), _id(++_id_generate) { }
      MGVertex(const MGVertex &old) : _atom(old._atom), _id(old._id) { }
      
      MGVertex& operator=(const MGVertex& other) {
        if (&other == this) return *this;
        _atom = other._atom;
        _id = other._id;
        return *this;
      }
      
      friend inline bool operator<(const MGVertex& lhs, const MGVertex& rhs) {
        return lhs._id < rhs._id;
      }
      
      friend inline bool operator==(const MGVertex& lhs, const MGVertex& rhs) {
        return !lhs._atom.expired() && lhs._atom.lock() == rhs._atom.lock();
      }
      
      Atom GetAtom() const { return _atom.lock(); }
      
      void SetAtom(Atom a) { _atom = a; }
    
    private:
      _Atom _atom;
      uid_ _id;
      static uid_ _id_generate;
    };
    
    class MGEdge {
    public:
      MGEdge() : _bond(), _id(++_id_generate) { }
      MGEdge(Bond b) : _bond(b), _id(++_id_generate) { }
      MGEdge(const MGEdge &old) : _bond(old._bond), _id(old._id) { }
      
      MGEdge& operator=(const MGEdge& other) {
        if (&other == this) return *this;
        _bond = other._bond;
        _id = other._id;
        return *this;
      }
      
      friend inline bool operator<(const MGEdge& lhs, const MGEdge& rhs) {
        return lhs._id < rhs._id;
      }
      
      friend inline bool operator==(const MGEdge& lhs, const MGEdge& rhs) {
        return !lhs._bond.expired() && lhs._bond.lock() == rhs._bond.lock();
      }
      
      Bond GetBond() const { return _bond.lock(); }
      void SetBond(Bond b) { _bond = b; }
      
    private:
      _Bond _bond;
      uid_ _id;
      static uid_ _id_generate;
    };
    
    class IXMolecularGraph : public GraphBase<MGVertex, MGEdge> {
    
    public:
      //! \cond
      IXMolecularGraph() = delete;
      //! \endcond
      
    private:
      //! \brief source molecule of the molecular graph.
      _Molecule _source;
    
    };
  }  // namespace graph
}  // namespace indigox

#endif /* INDIGOX_GRAPH_MOLECULAR_HPP */
