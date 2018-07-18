#ifndef INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP
#define INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/quad.hpp>
#include <indigox/utils/triple.hpp>

namespace indigox::test {
  struct IXAngle {
    indigox::IXAngle a;
    using AngAtmIter = indigox::IXAngle::AngleAtomIter;
    
    static Angle GetNewAngle() { return Angle(new indigox::IXAngle(Atom(), Atom(), Atom(), Molecule())); }
    static Angle GetNewAngle(Atom a, Atom b, Atom c, Molecule m) { return Angle(new indigox::IXAngle(a,b,c,m)); }
    
    // private wrapping members
    IXAngle() : a(Atom(), Atom(), Atom(), Molecule()) { }
    IXAngle(Atom i, Atom j, Atom k, Molecule m) : a(i,j,k,m) { }
    inline void Clear() { a.Clear(); }
    
    // public wrapping members
    inline stdx::triple<Atom,Atom,Atom> GetAtoms() { return a.GetAtoms(); }
    inline Molecule GetMolecule() { return a.GetMolecule(); }
    inline uid_ GetTag() { return a.GetTag(); }
    inline uid_ GetUniqueID() { return a.GetUniqueID(); }
    inline size_ NumAtoms() { return a.NumAtoms(); }
    inline void SetTag(uid_ t) { a.SetTag(t); }
    inline void SwapOrder() { a.SwapOrder(); }
    inline string_ ToString() { return a.ToString(); }
  };
  
  struct IXAtom {
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
    inline uid_ GetUniqueID() { return a.GetUniqueID(); }
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
    
  };
  
  struct IXBond {
    indigox::IXBond x;
  public:
    typedef indigox::IXBond::BondAtomIter BndAtmIter;
    
    static Bond GetNewBond() { return Bond(new indigox::IXBond(Atom(), Atom(), Molecule())); }
    static Bond GetNewBond(Atom a, Atom b, Molecule m) { return Bond(new indigox::IXBond(a,b,m)); }
    
    // private wrapping members
    IXBond() : x(Atom(), Atom(), Molecule()) { }
    IXBond(Atom a, Atom b, Molecule m) : x(a,b,m) { }
    inline void Clear() { x.Clear(); }
    
    // public wrapping members
    inline bool GetAromaticity() { return x.GetAromaticity(); }
    inline std::pair<Atom, Atom> GetAtoms() { return x.GetAtoms(); }
    inline Molecule GetMolecule() { return x.GetMolecule(); }
    inline BondOrder GetOrder() { return x.GetOrder(); }
    inline Atom GetSourceAtom() { return x.GetSourceAtom(); }
    inline BondStereo GetStereochemistry() { return x.GetStereochemistry(); }
    inline uid_ GetTag() { return x.GetTag(); }
    inline Atom GetTargetAtom() { return x.GetTargetAtom(); }
    inline uid_ GetUniqueID() { return x.GetUniqueID(); }
    inline size_ NumAtoms() { return x.NumAtoms(); }
    inline void SetAromaticity(bool a) { x.SetAromaticity(a); }
    inline void SetOrder(BondOrder a) { x.SetOrder(a); }
    inline void SetStereochemistry(BondStereo a) { x.SetStereochemistry(a); }
    inline void SetTag(uid_ tag) { x.SetTag(tag); }
    inline void SwapSourceTarget() { x.SwapSourceTarget(); }
    inline string_ ToString() { return x.ToString(); }
  };
  
  struct IXDihedral {
    indigox::IXDihedral x;
    using DhdAtmIter = indigox::IXDihedral::DihedAtomIter;
    
    static Dihedral GetNewDihedral() { return Dihedral(new indigox::IXDihedral(Atom(), Atom(), Atom(), Atom(), Molecule())); }
    static Dihedral GetNewDihedral(Atom a, Atom b, Atom c, Atom d, Molecule m) { return Dihedral(new indigox::IXDihedral(a,b,c,d,m)); }
    
    // private wrapping members
    IXDihedral() : x(Atom(), Atom(), Atom(), Atom(), Molecule()) { }
    IXDihedral(Atom a, Atom b, Atom c, Atom d, Molecule m) : x(a,b,c,d,m) { }
    inline void Clear() { x.Clear(); }
    
    // public wrapping members
    inline stdx::quad<Atom,Atom,Atom,Atom> GetAtoms() { return x.GetAtoms(); }
    inline Molecule GetMolecule() { return x.GetMolecule(); }
    inline uid_ GetTag() { return x.GetTag(); }
    inline uid_ GetUniqueId() { return x.GetUniqueID(); }
    inline size_ NumAtoms() { return x.NumAtoms(); }
    inline void SetTag(uid_ tag) { x.SetTag(tag); }
    inline void SwapOrder() { x.SwapOrder(); }
    inline string_ ToString() { return x.ToString(); }
  };
  
  
  struct IXMolecule {
    Molecule m; // needs to be a shared_ptr due to MolecularGraph ownership.
  public:
    typedef indigox::IXMolecule::MolAtomIter MolAtmIter;
    typedef indigox::IXMolecule::MolBondIter MolBndIter;
    typedef indigox::IXMolecule::MolAngleIter MolAngIter;
    typedef indigox::IXMolecule::MolDihedralIter MolDhdIter;
    typedef indigox::IXMolecule::Emergent Emergent;
    typedef std::bitset<static_cast<size_>(Emergent::NUM_EMERGENTS)> EmergeSet;
    // private wrapping members
    IXMolecule() : m(CreateMolecule()) { }
    inline void Init() { m->Init(); }
    inline Angle NewAngle(Atom a, Atom b, Atom c) { return m->NewAngle(a,b,c); }
    inline Dihedral NewDihedral(Atom a, Atom b, Atom c, Atom d) { return m->NewDihedral(a,b,c,d); }
    
    // public wrapping members
    inline Angle GetAngle(size_ pos) { return m->GetAngle(pos); }
    inline Angle GetAngle(Atom a, Atom b, Atom c) { return m->GetAngle(a,b,c); }
    inline Angle GetAngleID(uid_ t) { return m->GetAngleID(t); }
    inline std::pair<MolAngIter, MolAngIter> GetAngles() { return m->GetAngles(); }
    inline Angle GetAngleTag(uid_ tag) { return m->GetAngleTag(tag); }
    inline Atom GetAtom(size_ pos) { return m->GetAtom(pos); }
    inline Atom GetAtomID(uid_ id) { return m->GetAtomID(id); }
    inline std::pair<MolAtmIter, MolAtmIter> GetAtoms() { return m->GetAtoms(); }
    inline Atom GetAtomTag(uid_ tag) { return m->GetAtomTag(tag); }
    inline Bond GetBond(size_ pos) { return m->GetBond(pos); }
    inline Bond GetBond(Atom a, Atom b) { return m->GetBond(a,b); }
    inline Bond GetBondID(uid_ i) { return m->GetBondID(i); }
    inline std::pair<MolBndIter, MolBndIter> GetBonds() { return m->GetBonds(); }
    inline Bond GetBondTag(uid_ t) { return m->GetBondTag(t); }
    inline Dihedral GetDihedral(size_ p) { return m->GetDihedral(p); }
    inline Dihedral GetDihedral(Atom a, Atom b, Atom c, Atom d) { return m->GetDihedral(a,b,c,d); }
    inline Dihedral GetDihedralID(uid_ i) { return m->GetDihedralID(i); }
    inline std::pair<MolDhdIter, MolDhdIter> GetDihedrals() { return m->GetDihedrals(); }
    inline Dihedral GetDihedralTag(uid_ t) { return m->GetDihedralTag(t); }
    inline string_ GetFormula() { return m->GetFormula(); }
    inline graph::MolecularGraph GetGraph() { return m->GetGraph(); }
    inline int_ GetMolecularCharge() { return m->GetMolecularCharge(); }
    inline string_ GetName() { return m->GetName(); }
    inline uid_ GetUniqueID() { return m->GetUniqueID(); }
    inline bool HasAngle(Angle a) { return m->HasAngle(a); }
    inline bool HasAngle(Atom a, Atom b, Atom c) { return m->HasAngle(a,b,c); }
    inline bool HasAtom(Atom a) { return m->HasAtom(a); }
    inline bool HasBond(Bond b) { return m->HasBond(b); }
    inline bool HasBond(Atom a, Atom b) { return m->HasBond(a,b); }
    inline bool HasDihedral(Dihedral d) { return m->HasDihedral(d); }
    inline bool HasDihedral(Atom a, Atom b, Atom c, Atom d) { return m->HasDihedral(a,b,c,d); }
    inline Atom NewAtom() { return m->NewAtom(); }
    inline Atom NewAtom(Element e) { return m->NewAtom(e); }
    inline Atom NewAtom(string_ n) { return m->NewAtom(n); }
    inline Atom NewAtom(string_ n, Element e) { return m->NewAtom(n, e); }
    inline Bond NewBond(Atom a, Atom b) { return m->NewBond(a,b); }
    inline size_ NumAngles() { return m->NumAngles(); }
    inline size_ NumAtoms() { return m->NumAtoms(); }
    inline size_ NumBonds() { return m->NumBonds(); }
    inline size_ NumDihedrals() { return m->NumDihedrals(); }
    inline size_ PerceiveAngles() { return m->PerceiveAngles(); }
    inline size_ PerceiveDihedrals() { return m->PerceiveDihedrals(); }
    inline bool RemoveAtom(Atom a) { return m->RemoveAtom(a); }
    inline bool RemoveBond(Bond b) { return m->RemoveBond(b); }
    inline bool RemoveBond(Atom a, Atom b) { return m->RemoveBond(a,b); }
    inline void ReserveAtoms(size_ n) { m->ReserveAtoms(n); }
    inline void ReserveBonds(size_ n) { m->ReserveBonds(n); }
    inline void SetMolecularCharge(int_ q) { m->SetMolecularCharge(q); }
    inline void SetName(string_ n) { m->SetName(n); }
    inline void SetPropertyModified(MolProperty p) { m->SetPropertyModified(p); }
    
    // internals wrapping members
    // Emergent state
    inline EmergeSet GetEmergentState() { return m->_emerge; }
    inline void ResetEmergentState() { m->_emerge.reset(); }
    inline void ResetEmergentState(size_ p) { m->_emerge.reset(p); }
    inline void SetEmergentState(size_ p) { m->_emerge.set(p); }
    // Cached values
    inline string_ GetFormulaCache() { return m->_formula_cache; }
    // Container contains
    inline size_ AtomCapacity() { return m->_atms.capacity(); }
    inline size_ BondCapacity() { return m->_bnds.capacity(); }
    inline size_ CurrentAngleCount() { return m->_angs.size(); }
    inline size_ CurrentDihedralCount() { return m->_dhds.size(); }
  };
  
  
  struct IXAssignmentGraph {
    graph::IXAssignmentGraph g;
  public:
    using VertIter = graph::IXAssignmentGraph::VertIter;
    using NbrsIter = graph::IXAssignmentGraph::NbrsIter;
    // private wrapping members
    IXAssignmentGraph() : g(graph::MolecularGraph()) { }
    inline void AddEdges(graph::MGVertex s, graph::MGVertex t, graph::MGEdge e) { g.AddEdges(s,t,e); }
    inline graph::AGVertex AddVertex(graph::MGVertex v) { return g.AddVertex(v); }
    inline void DetermineAllNeighbours() { g.DetermineAllNeighbours(); }
    
    // public wrapping members
    IXAssignmentGraph(graph::MolecularGraph G) : g(G) { }
    inline size_ Degree(graph::AGVertex v) { return g.Degree(v); }
    inline std::pair<NbrsIter, NbrsIter> GetNeighbours(graph::AGVertex v) { return g.GetNeighbours(v); }
    inline graph::AGVertex GetVertex(graph::MGVertex v) { return g.GetVertex(v); }
    inline graph::AGVertex GetVertex(graph::MGEdge v) { return g.GetVertex(v); }
    inline std::pair<VertIter, VertIter> GetVertices() { return g.GetVertices(); }
    inline bool HasVertex(graph::AGVertex v) { return g.HasVertex(v); }
    inline bool HasVertex(graph::MGVertex v) { return g.HasVertex(v); }
    inline bool HasVertex(graph::MGEdge e) { return g.HasVertex(e); }
    inline bool IsConnected() { return g.IsConnected(); }
    inline size_ NumVertices() { return g.NumVertices(); }
  };
  
  struct IXMolecularGraph {
    graph::IXMolecularGraph g;
  public:
    typedef graph::IXMolecularGraph::EdgeIter EdgeIter;
    typedef graph::IXMolecularGraph::NbrsIter NbrsIter;
    typedef graph::IXMolecularGraph::VertIter VertIter;
    typedef graph::IXMolecularGraph::CompIter CompIter;
    
    // private wrapping members
    IXMolecularGraph() : g(Molecule()) { }
    IXMolecularGraph(Molecule m) : g(m) {}
    inline graph::MGEdge AddEdge(Bond bnd) { return g.AddEdge(bnd); }
    inline graph::MGVertex AddVertex(Atom atm) { return g.AddVertex(atm); }
    inline void Clear() { g.Clear(); }
    inline void RemoveEdge(graph::MGEdge e) { g.RemoveEdge(e); }
    inline void RemoveEdge(graph::MGVertex u, graph::MGVertex v) { g.RemoveEdge(u,v); }
    inline void RemoveVertex(graph::MGVertex v) { g.RemoveVertex(v); }
    
    // public wrapping methods
    inline size_ Degree(graph::MGVertex v) { return g.Degree(v); }
    inline std::pair<CompIter, CompIter> GetConnectedComponents() { return g.GetConnectedComponents(); }
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
    inline bool IsConnected() { return g.IsConnected(); }
    inline size_ NumConnectedComponents() { return g.NumConnectedComponents(); }
    inline size_ NumEdges() { return g.NumEdges(); }
    inline size_ NumVertices() { return g.NumVertices(); }
  };
  
  
}

#endif /* INDIGOX_TEST_CLASS_TEST_WRAPPERS_HPP */
