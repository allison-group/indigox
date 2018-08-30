#include <algorithm>

#include <EASTL/vector_set.h>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>

namespace indigox::graph {
  
  eastl::vector_set<Element> __con_elem = {
    GetPeriodicTable()->GetElement("H"),
    GetPeriodicTable()->GetElement("F"),
    GetPeriodicTable()->GetElement("Cl"),
    GetPeriodicTable()->GetElement("Br"),
    GetPeriodicTable()->GetElement("I")};

  test_suite_open("CondensedMolecularGraph");
  
// ============================================================================
// == IXCMGVertex CONSTRUCTION ================================================
// ============================================================================
  
  IXCMGVertex::IXCMGVertex(const MGVertex& v, const CondensedMolecularGraph& g)
  : _source(v), _graph(g) {
    using CS = ContractedSymmetry;
    MolecularGraph MG = g->GetSource();
    auto nbrs = MG->GetNeighbours(v);
    // Condensed vertices for the non-leaf addings
    if (std::distance(nbrs.first, nbrs.second) > 1) {
      for (; nbrs.first != nbrs.second; ++nbrs.first) {
        MGVertex u = *nbrs.first;
        MGEdge e = MG->GetEdge(u, v);
        if ((__con_elem.find(u->GetAtom()->GetElement()) != __con_elem.end())
            && (u->GetAtom()->GetFormalCharge() == 0)
            && (e->GetBond()->GetOrder() == BondOrder::SINGLE)) {
          Element el = u->GetAtom()->GetElement();
          CondensedVertex u_con;
          if (el == "H") u_con = std::make_pair(CS::Hydrogen, u);
          if (el == "F") u_con = std::make_pair(CS::Flourine, u);
          if (el == "Cl") u_con = std::make_pair(CS::Chlorine, u);
          if (el == "Br") u_con = std::make_pair(CS::Bromine, u);
          if (el == "I") u_con = std::make_pair(CS::Iodine, u);
          _con.insert(u_con);
        }
      }
    }
    
    // Determine Isomorphism Mask
    // Element number (7 bits)
    // Magnitude of formal charge (3 bits)
    // Is FC negative (1 bit)
    // Num condensed H (3 bits)
    // Num condensed F (3 bits)
    // Num condensed Cl (3 bits)
    // Num condensed Br (3 bits)
    // Num condensed I (3 bits)
    // Is in any cycle (1 bit)
    // Is in small cycle (<= 8)(1 bit)
    // Is R stereo (1 bit)
    // Is S stereo (1 bit)
    // Is aromatic (1 bit)
    VertexIsoMask atm_num, fc_mag, h, f, cl, br, i;
    atm_num.from_uint64(v->GetAtom()->GetElement()->GetAtomicNumber());
    fc_mag.from_uint64(abs(v->GetAtom()->GetFormalCharge())); fc_mag <<= 7;
    h.from_uint32(NumContracted(CS::Hydrogen)); h <<= 11;
    f.from_uint32(NumContracted(CS::Flourine)); f <<= 14;
    cl.from_uint32(NumContracted(CS::Chlorine)); cl <<= 17;
    br.from_uint32(NumContracted(CS::Bromine)); br <<= 20;
    i.from_uint32(NumContracted(CS::Iodine)); i <<= 23;
    _iso_mask = atm_num | fc_mag | h | f | cl | br | i;
    if (v->GetAtom()->GetFormalCharge() < 0) _iso_mask.set(10);
//    if (v->IsCyclic()) _iso_mask.set(26);
//    if (v->IsCyclic(8)) _iso_mask.set(27);
    if (v->GetAtom()->GetStereochemistry() == AtomStereo::R) _iso_mask.set(28);
    if (v->GetAtom()->GetStereochemistry() == AtomStereo::S) _iso_mask.set(29);
    if (v->GetAtom()->GetAromaticity()) _iso_mask.set(30);
  }
  
  size_ IXCMGVertex::NumContracted(ContractedSymmetry sym) const {
    auto counter = [&sym](const CondensedVertex& v) { return v.first == sym; };
    return std::count_if(_con.begin(), _con.end(), counter);
  }
  
  bool IXCMGVertex::IsContractedHere(const MGVertex& v) const {
    auto finder = [&v](const CondensedVertex& u) { return u.second == v; };
    return std::find_if(_con.begin(), _con.end(), finder) != _con.end();
  }
  
  
// ============================================================================
// == IXCMGEdge CONSTRUCTION ==================================================
// ============================================================================
  
  IXCMGEdge::IXCMGEdge(const MGEdge& e, const CondensedMolecularGraph& g)
  : _source(e), _graph(g) {
    // Determine Isomorphism Mask
    // Order (3 bits)
    // Is E stereo (1 bit)
    // Is Z stereo (1 bit)
    // Is in any cycle (1 bit)
    // Is in small cycle (<= 8)(1 bit)
    // Is aromatic (1 bit)
    _iso_mask.from_uint32(static_cast<uint32_t>(e->GetBond()->GetOrder()));
    if (e->GetBond()->GetStereochemistry() == BondStereo::E) _iso_mask.set(3);
    if (e->GetBond()->GetStereochemistry() == BondStereo::Z) _iso_mask.set(4);
//    if (e->IsCyclic()) _iso_mask.set(5);
//    if (e->IsCyclic(8)) _iso_mask.set(6);
    if (e->GetBond()->GetAromaticity()) _iso_mask.set(7);
  }
  
  
// ============================================================================
// == IXCondensedMolecularGraph CONSTRUCTION ==================================
// ============================================================================
  
  IXCondensedMolecularGraph::IXCondensedMolecularGraph(const MolecularGraph& g)
  : _source(g->CreateSnapshot()), _g() { }
  
  CondensedMolecularGraph CondenseMolecularGraph(const MolecularGraph& G) {
    CondensedMolecularGraph G2 = std::make_shared<IXCondensedMolecularGraph>(G);
    MolecularGraph Gsource = G2->_source;
    auto vertices = Gsource->GetVertices();
    for (; vertices.first != vertices.second; ++vertices.first) {
      MGVertex v = *vertices.first;
      // Add vertex for all non-leaf vertices
      if (Gsource->Degree(v) > 1) { G2->AddVertex(v); continue; }
      Element e = v->GetAtom()->GetElement();
      // Add vertex if leaf vertex is not in __con_elem
      if (__con_elem.find(e) == __con_elem.end()) { G2->AddVertex(v); continue; }
      // Add vertex if has formal charge
      if (v->GetAtom()->GetFormalCharge() != 0) { G2->AddVertex(v); continue; }
      // Add vertex if bond is not single
      MGVertex u = *(Gsource->GetNeighbours(v).first);
      BondOrder order = Gsource->GetEdge(v, u)->GetBond()->GetOrder();
      if (order == BondOrder::UNDEFINED) {
        throw std::runtime_error("Attempting to condense a molecule with undefined bond orders");
      }
      if (order != BondOrder::SINGLE) {
        G2->AddVertex(v); continue; }
      // Add vertex if parent is also a leaf
      if (Gsource->Degree(u) == 1) G2->AddVertex(v);
    }
    
    auto edges = Gsource->GetEdges();
    for (; edges.first != edges.second; ++edges.first) {
      MGEdge e = *edges.first;
      MGVertex u, v;
      std::tie(u,v) = Gsource->GetVertices(e);
      if (G2->HasVertex(u) && G2->HasVertex(v)) G2->AddEdge(e);
    }
    
    return G2;
  }
  
  CMGEdge IXCondensedMolecularGraph::AddEdge(const MGEdge &e) {
    CMGEdge E = std::make_shared<IXCMGEdge>(e, shared_from_this());
    CMGVertex u = GetVertex(_source->GetSource(e));
    CMGVertex v = GetVertex(_source->GetTarget(e));
    _g.AddEdge(u.get(), v.get(), E.get());
    _emap.emplace(e, E);
    _e.emplace_back(E);
    _n[u].emplace_back(v);
    _n[v].emplace_back(u);
    return E;
  }
  
  CMGVertex IXCondensedMolecularGraph::AddVertex(const MGVertex &v) {
    CMGVertex V = std::make_shared<IXCMGVertex>(v, shared_from_this());
    _g.AddVertex(V.get());
    _vmap.emplace(v, V);
    _v.emplace_back(V);
    _n.emplace(V, NbrsContain::mapped_type());
    return V;
  }
  
  test_suite_close();
}
