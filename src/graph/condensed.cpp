#include <algorithm>

#include <EASTL/vector_set.h>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>

namespace indigox::graph {
  
  eastl::vector_set<Element> __con_elem = {
    GetPeriodicTable().GetElement("H"),
    GetPeriodicTable().GetElement("F"),
    GetPeriodicTable().GetElement("Cl"),
    GetPeriodicTable().GetElement("Br"),
    GetPeriodicTable().GetElement("I")};

  test_suite_open("CondensedMolecularGraph");
  
  struct CMGVertex::CMGVertexData {
    MGVertex source;
    wCondensedMolecularGraph graph;
    std::vector<CondensedVertex> condensed;
    VertexIsoMask mask;
    
    CMGVertexData() = default;
    CMGVertexData(const MGVertex& v, CondensedMolecularGraph& g)
    : source(v), graph(g.weak_from_this()), mask(0) { }
    
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("source", source),
              INDIGOX_SERIAL_NVP("graph", graph),
              INDIGOX_SERIAL_NVP("condensed_verts", condensed),
              INDIGOX_SERIAL_NVP("mask", mask));
    }
  };
  
  struct CMGEdge::CMGEdgeData {
    MGEdge source;
    wCondensedMolecularGraph graph;
    EdgeIsoMask mask;
    
    CMGEdgeData() = default;
    CMGEdgeData(const MGEdge& e, CondensedMolecularGraph& g)
    : source(e), graph(g.weak_from_this()), mask(0) { }
    
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("source", source),
              INDIGOX_SERIAL_NVP("graph", graph),
              INDIGOX_SERIAL_NVP("mask", mask));
    }
  };
  
// ============================================================================
// == CMGVertex CONSTRUCTION ==================================================
// ============================================================================
  
  CMGVertex::CMGVertex() : _dat(nullptr) { }
  CMGVertex::CMGVertex(const CMGVertex& v) : _dat(v._dat) { }
  CMGVertex::CMGVertex(CMGVertex&& v) noexcept : _dat(std::move(v._dat)) { }
  CMGVertex& CMGVertex::operator=(const CMGVertex &v) {
    if (&v != this) _dat = v._dat;
    return *this;
  }
  CMGVertex& CMGVertex::operator=(CMGVertex &&v) {
    _dat = std::move(v._dat);
    return *this;
  }
  
  CMGVertex::ContractedSymmetry __get_symmetry(size_t atomic_number,
                                               const MGVertex&) {
    using CS = CMGVertex::ContractedSymmetry;
    switch (atomic_number) {
      case 1: return CS::Hydrogen;
      case 9: return CS::Fluorine;
      case 17: return CS::Chlorine;
      case 35: return CS::Bromine;
      case 53: return CS::Iodine;
      default: throw std::runtime_error("Cannot determine contracted symmetry");
    }
  }
  
  CMGVertex::CMGVertex(const MGVertex& v, CondensedMolecularGraph& g)
  : _dat(std::make_shared<CMGVertexData>(v, g)) {
    using CS = ContractedSymmetry;
    MolecularGraph& MG = g.GetMolecularGraph();
    const MolecularGraph::VertContain& nbrs = MG.GetNeighbours(v);
    // Condensed vertices for the non-leaf addings
    if (nbrs.size() > 1) {
      for (const MGVertex& u : nbrs) {
        MGEdge e = MG.GetEdge(u, v);
        Atom& atm = u.GetAtom();
        Bond& bnd = e.GetBond();
        if ((__con_elem.find(atm.GetElement()) != __con_elem.end())
            && (atm.GetFormalCharge() == 0)
            && (bnd.GetOrder() == BondOrder::SINGLE)) {
          size_t el = atm.GetElement().GetAtomicNumber();
          CondensedVertex u_con = std::make_pair(__get_symmetry(el, u), u);
          _dat->condensed.emplace_back(u_con);
        }
      }
      std::sort(_dat->condensed.begin(), _dat->condensed.end());
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
    Atom& atm = v.GetAtom();
    VertexIsoMask atm_num, fc_mag, h, f, cl, br, i, mask, degree;
    atm_num.from_uint64(atm.GetElement().GetAtomicNumber());
    fc_mag.from_uint64(abs(atm.GetFormalCharge())); fc_mag <<= 7;
    h.from_uint32(NumContracted(CS::Hydrogen) + atm.GetImplicitCount()); h <<= 11;
    f.from_uint32(NumContracted(CS::Fluorine)); f <<= 14;
    cl.from_uint32(NumContracted(CS::Chlorine)); cl <<= 17;
    br.from_uint32(NumContracted(CS::Bromine)); br <<= 20;
    i.from_uint32(NumContracted(CS::Iodine)); i <<= 23;
    degree.from_uint32(atm.NumBonds()); degree <<= 31;
    mask = atm_num | fc_mag | h | f | cl | br | i | degree;
    if (atm.GetFormalCharge() < 0) mask.set(10);
    if (MG.IsCyclic(v, 8)) mask.set(26);
//    if (v->IsCyclic(8)) _iso_mask.set(27);
    if (atm.GetStereochemistry() == AtomStereo::R) mask.set(28);
    if (atm.GetStereochemistry() == AtomStereo::S) mask.set(29);
    if (atm.GetAromaticity()) mask.set(30);
    
    _dat->mask = mask;
  }
  
// ============================================================================
// == CMGVertex Serialisation =================================================
// ============================================================================
  template <typename Archive>
  void CMGVertex::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", _dat));
  }
  INDIGOX_SERIALISE(CMGVertex);
  
// ============================================================================
// == CMGVertex Methods =======================================================
// ============================================================================
  MGVertex CMGVertex::GetSource() const { return _dat->source; }
  
  CondensedMolecularGraph& CMGVertex::GetGraph() const {
    return *_dat->graph.lock();
  }
  
  size_t CMGVertex::NumContracted() const { return _dat->condensed.size(); }
  
  size_t CMGVertex::NumContracted(ContractedSymmetry s) const {
    return std::accumulate(_dat->condensed.begin(), _dat->condensed.end(), 0,
                           [&s](size_t a, const CondensedVertex& v) -> size_t {
                             return s == v.first ? ++a : a; });
  }
  
  VertexIsoMask& CMGVertex::GetIsomorphismMask() const { return _dat->mask; }
  
  bool CMGVertex::IsContractedHere(const MGVertex &v) const {
    size_t z = v.GetAtom().GetElement().GetAtomicNumber();
    CondensedVertex cv = std::make_pair(__get_symmetry(z, v), v);
    return std::find(_dat->condensed.begin(), _dat->condensed.end(), cv)
            != _dat->condensed.end();
  }
  
  std::vector<MGVertex> CMGVertex::GetContractedVertices() const {
    std::vector<MGVertex> verts;
    verts.reserve(NumContracted());
    for (auto& cv : _dat->condensed) verts.emplace_back(cv.second);
    return verts;
  }
  
  const std::vector<CMGVertex::CondensedVertex>&
  CMGVertex::GetCondensedVertices() const {
    return _dat->condensed;
  }
  
  std::ostream& operator<<(std::ostream& os, const CMGVertex& v) {
    if (v) os << "CMGVertex(" << v.GetSource().GetAtom().GetTag() << ")";
    return os;
  }
  
// ============================================================================
// == CMGEdge CONSTRUCTION ====================================================
// ============================================================================
  
  CMGEdge::CMGEdge() : _dat(nullptr) { }
  CMGEdge::CMGEdge(const CMGEdge& e) : _dat(e._dat) { }
  CMGEdge::CMGEdge(CMGEdge&& e) noexcept : _dat(e._dat) { }
  CMGEdge& CMGEdge::operator=(const CMGEdge &e) {
    if (&e != this) _dat = e._dat;
    return *this;
  }
  CMGEdge& CMGEdge::operator=(CMGEdge &&e) {
    _dat = std::move(e._dat);
    return *this;
  }
  
  CMGEdge::CMGEdge(const MGEdge& e, CondensedMolecularGraph& g)
  : _dat(std::make_shared<CMGEdgeData>(e, g)) {
    EdgeIsoMask mask, degree_small, degree_large;
    // Determine Isomorphism Mask
    // Order (3 bits)
    // Is E stereo (1 bit)
    // Is Z stereo (1 bit)
    // Is in any cycle (1 bit)
    // Is in small cycle (<= 8)(1 bit)
    // Is aromatic (1 bit)
    Bond& bnd = e.GetBond();
    Atom& a = bnd.GetSourceAtom();
    Atom& b = bnd.GetTargetAtom();
    degree_small.from_uint32(a.NumBonds()); degree_small <<= 8;
    degree_large.from_uint32(b.NumBonds()); degree_large <<= 11;
    if (a.NumBonds() > b.NumBonds()) std::swap(degree_small, degree_large);
    mask.from_uint32(static_cast<uint32_t>(bnd.GetOrder()));
    mask |= degree_small | degree_large;
    if (bnd.GetStereochemistry() == BondStereo::E) mask.set(3);
    if (bnd.GetStereochemistry() == BondStereo::Z) mask.set(4);
    if (g.GetMolecularGraph().IsCyclic(e, 8)) mask.set(5);
//    if (e->IsCyclic(8)) _iso_mask.set(6);
    if (bnd.GetAromaticity()) mask.set(7);
    _dat->mask = mask;
  }
  
  
// ============================================================================
// == CMGEdge Serialisation ===================================================
// ============================================================================
  template <typename Archive>
  void CMGEdge::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", _dat));
  }
  INDIGOX_SERIALISE(CMGEdge);
  
// ============================================================================
// == CMGEdge Methods =========================================================
// ============================================================================
  
  MGEdge CMGEdge::GetSource() const { return _dat->source; }
  CondensedMolecularGraph& CMGEdge::GetGraph() const { return *_dat->graph.lock(); }
  EdgeIsoMask& CMGEdge::GetIsomorphismMask() const { return _dat->mask; }

// ============================================================================
// == CondensedMolecularGraph Serialisation ===================================
// ============================================================================
  
  template <typename Archive>
  void CondensedMolecularGraph::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("graph", cereal::base_class<graph_type>(this)),
            INDIGOX_SERIAL_NVP("verts", _vmap),
            INDIGOX_SERIAL_NVP("edges", _emap),
            INDIGOX_SERIAL_NVP("source", _source),
            INDIGOX_SERIAL_NVP("supergraph", _subg));
  }
  template <typename Archive>
  void CondensedMolecularGraph::load(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("graph", cereal::base_class<graph_type>(this)),
            INDIGOX_SERIAL_NVP("verts", _vmap),
            INDIGOX_SERIAL_NVP("edges", _emap),
            INDIGOX_SERIAL_NVP("source", _source),
            INDIGOX_SERIAL_NVP("supergraph", _subg));
  }
  INDIGOX_SERIALISE_SPLIT(CondensedMolecularGraph);

// ============================================================================
// == CondensedMolecularGraph Modification ====================================
// ============================================================================
  
  CMGEdge CondensedMolecularGraph::AddEdge(const MGEdge &e) {
    MolecularGraph& MG = e.GetGraph();
    CMGVertex u = GetVertex(MG.GetSourceVertex(e));
    CMGVertex v = GetVertex(MG.GetTargetVertex(e));
    CMGEdge e_local = CMGEdge(e, *this);
    _emap.emplace(e, e_local);
    graph_type::AddEdge(u, v, e_local);
    return e_local;
  }
  
  CMGVertex CondensedMolecularGraph::AddVertex(const MGVertex &v) {
    CMGVertex v_locl = CMGVertex(v, *this);
    _vmap.emplace(v, v_locl);
    graph_type::AddVertex(v_locl);
    return v_locl;
  }
  
  void CondensedMolecularGraph::Clear() {
//    _vmap.clear();
//    _emap.clear();
//    _source.reset();
//    graph_type::Clear();
  }
  
// ============================================================================
// == CondensedMolecularGraph CONSTRUCTION ====================================
// ============================================================================
  
  CondensedMolecularGraph::CondensedMolecularGraph()
  : _source(wMolecularGraph()), _subg() { }
  
  CondensedMolecularGraph::CondensedMolecularGraph(MolecularGraph& g)
  : _source(g.weak_from_this()), _subg() { }
  
  sCondensedMolecularGraph Condense(MolecularGraph& MG) {
    using CMG = sCondensedMolecularGraph;
    CMG CG = CMG(new CondensedMolecularGraph(MG));
    
    for (const MGVertex& v : MG.GetVertices()) {
      // Add vertex for all non-leaf vertices
      if (MG.Degree(v) > 1) { CG->AddVertex(v); continue; }
      Atom& atm = v.GetAtom();
      Element e = atm.GetElement();
      // Add vertex if leaf vertex is not in __con_elem
      if (__con_elem.find(e) == __con_elem.end()) { CG->AddVertex(v); continue; }
      // Add vertex if has formal charge
      if (atm.GetFormalCharge() != 0) { CG->AddVertex(v); continue; }
      // Add vertex if bond is not single
      MGVertex u = MG.GetNeighbours(v).front();
      BondOrder order = MG.GetEdge(v, u).GetBond().GetOrder();
      if (order != BondOrder::SINGLE) { CG->AddVertex(v); continue; }
      // Add vertex if parent is also a leaf
      if (MG.Degree(u) == 1) CG->AddVertex(v);
    }
    
    for (const MGEdge& e : MG.GetEdges()) {
      MGVertex u, v;
      std::tie(u,v) = MG.GetVertices(e);
      if (CG->HasVertex(u) && CG->HasVertex(v)) CG->AddEdge(e);
    }
    return CG;
  }
  
// ============================================================================
// == CondensedMolecularGraph Getting and Checking ============================
// ============================================================================
  MolecularGraph& CondensedMolecularGraph::GetMolecularGraph() {
    if (IsSubgraph())
      throw std::runtime_error("Subgraphs do not relate to molecular graphs");
    return *_source.lock();
  }
  
  CondensedMolecularGraph& CondensedMolecularGraph::GetSuperGraph() {
    if (!IsSubgraph())
      throw std::runtime_error("Cannot get supergraph of non subgraph");
    return *_subg;
  }
  
  CMGEdge CondensedMolecularGraph::GetEdge(const MGEdge &e) const {
    return _emap.at(e); }
  
  CMGVertex CondensedMolecularGraph::GetVertex(const MGVertex &v) const {
    return _vmap.at(v); }
  
  CMGVertex CondensedMolecularGraph::GetCondensedVertex(const MGVertex &v) const {
    return *std::find_if(_v.begin(), _v.end(), [&v](const CMGVertex& u) {
      return u.IsContractedHere(v); });
  }
  
  bool CondensedMolecularGraph::HasVertex(const MGVertex &v) const {
    return _vmap.find(v) != _vmap.end(); }
  
  bool CondensedMolecularGraph::HasEdge(const MGEdge &e) const {
    return _emap.find(e) != _emap.end(); }
  
  bool CondensedMolecularGraph::HasCondensedVertex(const MGVertex &v) const {
    return std::find_if(_v.begin(), _v.end(), [&v](const CMGVertex& u) {
      return u.IsContractedHere(v); }) != _v.end(); }
  
// ============================================================================
// == CondensedMolecularGraph Subgraph generation =============================
// ============================================================================
  sCondensedMolecularGraph
  CondensedMolecularGraph::Subgraph(std::vector<CMGVertex> &verts) {
    sCondensedMolecularGraph G = std::make_shared<CondensedMolecularGraph>();
    G->_subg = shared_from_this();
    for (const CMGVertex& v : verts) {
      if (!HasVertex(v)) throw std::runtime_error("Non-member vertex");
      MGVertex source = v.GetSource();
      G->_vmap.emplace(source, v);
      G->graph_type::AddVertex(v);
    }
    
    for (const CMGEdge& e : GetEdges()) {
      CMGVertex u = GetSourceVertex(e);
      CMGVertex v = GetTargetVertex(e);
      if (!G->HasVertex(u) || !G->HasVertex(v)) continue;
      G->_emap.emplace(e.GetSource(), e);
      G->graph_type::AddEdge(u, v, e);
    }
    return G;
  }
  
  sCondensedMolecularGraph
  CondensedMolecularGraph::Subgraph(std::vector<CMGVertex> &verts,
                                    std::vector<CMGEdge> &edges) {
    sCondensedMolecularGraph G = std::make_shared<CondensedMolecularGraph>();
    G->_subg = shared_from_this();
    for (const CMGVertex& v : verts) {
      if (!HasVertex(v)) throw std::runtime_error("Non-member vertex");
      MGVertex source = v.GetSource();
      G->_vmap.emplace(source, v);
      G->graph_type::AddVertex(v);
    }
    
    for (const CMGEdge& e : edges) {
      if (!HasEdge(e)) throw std::runtime_error("Non-member edge");
      CMGVertex u = GetSourceVertex(e);
      CMGVertex v = GetTargetVertex(e);
      if (!G->HasVertex(u) || !G->HasVertex(v)) continue;
      G->_emap.emplace(e.GetSource(), e);
      G->graph_type::AddEdge(u, v, e);
    }
    return G;
  }
  
  test_suite_close();
}
