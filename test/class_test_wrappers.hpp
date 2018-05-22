#ifndef INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP
#define INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
//#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/numerics.hpp>

// Dummy classes
namespace indigox {
  class IXAngle { };
  class IXDihedral { };
  class IXMolecule { };
}

namespace indigox::test {
  class IXAtom {
    indigox::IXAtom a;
  public:
    typedef indigox::IXAtom::AtomAngleIter AtmAngIter;
    typedef indigox::IXAtom::AtomBondIter AtmBndIter;
    typedef indigox::IXAtom::AtomDihedralIter AtmDhdIter;
    
    static Atom GetNewAtom() { return Atom(new indigox::IXAtom(Molecule())); }
    static Atom GetNewAtom(Molecule mol) { return Atom(new indigox::IXAtom(mol)); }
    
    // private wrapping members
    IXAtom() : a(Molecule()) { }
    IXAtom(Molecule m) : a(m) { }
    inline void AddAngle(Angle i) { a.AddAngle(i); }
    inline void AddBond(Bond i) { a.AddBond(i); }
    inline void AddDihedral(Dihedral i) { a.AddDihedral(i); }
    inline void Clear() { a.Clear(); }
    inline void RemoveAngle(Angle i) { a.RemoveAngle(i); }
    inline void RemoveBond(Bond i) { a.RemoveBond(i); }
    inline void RemoveDihedral(Dihedral i) { a.RemoveDihedral(i); }
    
    // public wrapping members
    inline uint_ AddImplicitHydrogen() { return a.AddImplicitHydrogen(); }
    inline std::pair<AtmAngIter, AtmAngIter> GetAngleIters() { return a.GetAngleIters(); }
    inline bool GetAromaticity() { return a.GetAromaticity(); }
    inline std::pair<AtmBndIter, AtmBndIter> GetBondIters() { return a.GetBondIters(); }
    inline std::pair<AtmDhdIter, AtmDhdIter> GetDihedralIters() { return a.GetDihedralIters(); }
    inline Element GetElement() { return a.GetElement(); }
    inline int_ GetFormalCharge() { return a.GetFormalCharge(); }
    inline uint_ GetImplicitCount() { return a.GetImplicitCount(); }
    inline Molecule GetMolecule() { return a.GetMolecule(); }
    inline string_ GetName() { return a.GetName(); }
    inline float_ GetPartialCharge() { return a.GetPartialCharge(); }
    inline AtomStereo GetStereochemistry() { return a.GetStereochemistry(); }
    inline uint_ GetTag() { return a.GetTag(); }
    inline Vec3 GetVector() { return a.GetVector(); }
    inline float_ GetX() { return a.GetX(); }
    inline float_ GetY() { return a.GetY(); }
    inline float_ GetZ() { return a.GetZ(); }
    inline size_ NumAngles() { return a.NumAngles(); }
    inline size_ NumBonds() { return a.NumBonds(); }
    inline size_ NumDihedrals() { return a.NumDihedrals(); }
    inline uint_ RemoveImplicitHydrogen() { return a.RemoveImplicitHydrogen(); }
    inline void SetAromaticity(bool i) { a.SetAromaticity(i); }
    inline void SetElement(Element i) { a.SetElement(i); }
    inline void SetElement(string_ i) { a.SetElement(i); }
    inline void SetElement(uint_ i) { a.SetElement(i); }
    inline void SetFormalCharge(int_ i) { a.SetFormalCharge(i); }
    inline void SetImplicitCount(uint_ i) { a.SetImplicitCount(i); }
    inline void SetName(string_ i) { a.SetName(i); }
    inline void SetPartialCharge(float_ i) { a.SetPartialCharge(i); }
    inline void SetPosition(float_ x, float_ y, float_ z) { a.SetPosition(x,y,z); }
    inline void SetStereochemistry(AtomStereo i) { a.SetStereochemistry(i); }
    inline void SetTag(uint_ i) { a.SetTag(i); }
    inline void SetX(float_ i) { a.SetX(i); }
    inline void SetY(float_ i) { a.SetY(i); }
    inline void SetZ(float_ i) { a.SetZ(i); }
    inline string_ ToString() { return a.ToString(); }
    inline uid_ GetUniqueID() { return a.GetUniqueID(); }
    
  };
  
  class IXBond {
    indigox::IXBond x;
  public:
    typedef indigox::IXBond::BondAtomIter BndAtmIter;
    typedef indigox::IXBond::BondAngleIter BndAngIter;
    typedef indigox::IXBond::BondDihedralIter BndDhdIter;
    
    // private wrapping members
    IXBond() : x(Atom(), Atom(), Molecule()) { }
    IXBond(Atom a, Atom b, Molecule m) : x(a,b,m) { }
    inline void AddAngle(Angle a) { x.AddAngle(a); }
    inline void AddDihedral(Dihedral a) { x.AddDihedral(a); }
    inline void Clear() { x.Clear(); }
    inline void RemoveAngle(Angle a) { x.RemoveAngle(a); }
    inline void RemoveDihedral(Dihedral a) { x.RemoveDihedral(a); }
    
    // public wrapping members
    inline std::pair<BndAngIter, BndAngIter> GetAngleIters() { return x.GetAngleIters(); }
    inline bool GetAromaticity() { return x.GetAromaticity(); }
    inline std::pair<BndAtmIter, BndAtmIter> GetAtomIters() { return x.GetAtomIters(); }
    inline std::pair<BndDhdIter, BndDhdIter> GetDihedralIters() { return x.GetDihedralIters(); }
    inline Molecule GetMolecule() { return x.GetMolecule(); }
    inline BondOrder GetOrder() { return x.GetOrder(); }
    inline Atom GetSourceAtom() { return x.GetSourceAtom(); }
    inline BondStereo GetStereochemistry() { return x.GetStereochemistry(); }
    inline uid_ GetTag() { return x.GetTag(); }
    inline Atom GetTargetAtom() { return x.GetTargetAtom(); }
    inline size_ NumAngles() { return x.NumAngles(); }
    inline size_ NumAtoms() { return x.NumAtoms(); }
    inline size_ NumDihedrals() { return x.NumDihedrals(); }
    inline void SetAromaticity(bool a) { x.SetAromaticity(a); }
    inline void SetOrder(BondOrder a) { x.SetOrder(a); }
    inline void SetStereochemistry(BondStereo a) { x.SetStereochemistry(a); }
    inline void SetTag(uid_ tag) { x.SetTag(tag); }
    inline void SwapSourceTarget() { x.SwapSourceTarget(); }
    inline string_ ToString() { return x.ToString(); }
    inline uid_ GetUniqueID() { return x.GetUniqueID(); }
  };
  
  class IXMolecularGraph {
    graph::IXMolecularGraph g;
  public:
    typedef graph::IXMolecularGraph::EdgeIter EdgeIter;
    typedef graph::IXMolecularGraph::NbrsIter NbrsIter;
    typedef graph::IXMolecularGraph::VertIter VertIter;
    IXMolecularGraph() = delete;
    IXMolecularGraph(Molecule m) : g(m) {}
    inline size_ Degree(graph::MGVertex v) { return g.Degree(v); }
    inline graph::MGEdge GetEdge(graph::MGVertex u, graph::MGVertex v) { return g.GetEdge(u, v); }
    inline graph::MGEdge GetEdge(Bond b) { return g.GetEdge(b); }
    inline std::pair<EdgeIter, EdgeIter> GetEdges() { return g.GetEdges(); }
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(graph::MGVertex v) { return g.GetNeighbours(v); }
    inline graph::MGVertex GetSource(graph::MGEdge e) { return g.GetSource(e); }
    inline graph::MGVertex GetTarget(graph::MGEdge e) { return g.GetTarget(e); }
    inline graph::MGVertex GetVertex(Atom a) { return g.GetVertex(a); }
    inline std::pair<graph::MGVertex, graph::MGVertex> GetVertices(graph::MGEdge e) { return g.GetVertices(e); }
    inline std::pair<VertIter, VertIter> GetVertices() { return g.GetVertices(); }
    inline bool HasEdge(Bond b) { return g.HasEdge(b); }
    inline bool HasEdge(graph::MGEdge e) { return g.HasEdge(e); }
    inline bool HasEdge(graph::MGVertex u, graph::MGVertex v) { return g.HasEdge(u,v); }
    inline bool HasVertex(Atom v) { return g.HasVertex(v); }
    inline bool HasVertex(graph::MGVertex v) { return g.HasVertex(v); }
    inline size_ NumEdges() { return g.NumEdges(); }
    inline size_ NumVertices() { return g.NumVertices(); }
    //
    inline graph::MGEdge AddEdge(Bond bnd) { return g.AddEdge(bnd); }
    inline graph::MGVertex AddVertex(Atom atm) { return g.AddVertex(atm); }
    inline void Clear() { g.Clear(); }
    inline void RemoveEdge(graph::MGEdge e) { g.RemoveEdge(e); }
    inline void RemoveEdge(graph::MGVertex u, graph::MGVertex v) { g.RemoveEdge(u,v); }
    inline void RemoveVertex(graph::MGVertex v) { g.RemoveVertex(v); }
  };
}

#endif /* INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP */
