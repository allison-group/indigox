#include <algorithm>
#include <limits>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/molecular_test.hpp>
#include <indigox/test/atom_test.hpp>
#include <indigox/test/bond_test.hpp>

namespace indigox::graph {
  test_suite_open("MolecularGraph");
  
  struct MGVertex::MGVertexData {
    wAtom atom;
    wMolecularGraph graph;
    
    MGVertexData() = default;
    MGVertexData(wAtom a, wMolecularGraph g) : atom(a), graph(g) { }
    
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("atom", atom),
              INDIGOX_SERIAL_NVP("graph", graph));
    }
  };
  
  struct MGEdge::MGEdgeData {
    wBond bond;
    wMolecularGraph graph;
    
    MGEdgeData() = default;
    MGEdgeData(wBond b, wMolecularGraph g) : bond(b), graph(g) { }
    
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("bond", bond),
              INDIGOX_SERIAL_NVP("graph", graph));
    }
  };
  
// ============================================================================
// == SERIALISATION ===========================================================
// ============================================================================
  
  template <typename Archive>
  void MGVertex::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", _dat));
  }
  INDIGOX_SERIALISE(MGVertex);
  
/*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMGVertex serialisation", T, ixmolv_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestMolecularVertex saved(test::CreateGenericTestAtom().imp,
                                    test::CreateGenericTestMolecularGraph().imp);
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp, saved.GetAtom(), saved.GetGraph()));
    }
    
    Atom atom_loaded;
    graph::MGVertex loaded;
    graph::MolecularGraph graph_loaded;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, atom_loaded, graph_loaded));
    }
    
    check_eq(loaded->GetAtom(), atom_loaded);
    check_eq(loaded->GetGraph(), graph_loaded);
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixmolv_serial, ixserial<IXMGVertex>);
 */
  
  template <typename Archive>
  void MGEdge::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", _dat));
  }
  INDIGOX_SERIALISE(MGEdge);
  
/*  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMGEdge serialisation", T, ixmole_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestMolecularEdge saved(test::CreateGenericTestBond().imp,
                                  test::CreateGenericTestMolecularGraph().imp);
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp, saved.GetBond(), saved.GetGraph()));
    }
    
    Bond bond_loaded;
    MGEdge loaded;
    MolecularGraph loaded_graph;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, bond_loaded, loaded_graph));
    }
    
    check_eq(loaded->GetBond(), bond_loaded);
    check_eq(loaded->GetGraph(), loaded_graph);
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixmole_serial, ixserial<IXMGEdge>);
 */
  
  template <typename Archive>
  void MolecularGraph::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("graph", cereal::base_class<graph_type>(this)),
            INDIGOX_SERIAL_NVP("atom_map", _at2v),
            INDIGOX_SERIAL_NVP("bond_map", _bn2e),
            INDIGOX_SERIAL_NVP("molecule", _mol),
            INDIGOX_SERIAL_NVP("condensed_graph", _cond));
  }
  
  template <typename Archive>
  void MolecularGraph::load(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("graph", cereal::base_class<graph_type>(this)),
            INDIGOX_SERIAL_NVP("atom_map", _at2v),
            INDIGOX_SERIAL_NVP("bond_map", _bn2e),
            INDIGOX_SERIAL_NVP("molecule", _mol),
            INDIGOX_SERIAL_NVP("condensed_graph", _cond));
  }
  
  INDIGOX_SERIALISE_SPLIT(MolecularGraph);
  
 /* DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMolecularGraph serialisation", T, ixmolg_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    test::MolecularGraphTestFixture fixture;
    
    MolecularGraph saved = fixture.benzene.imp;
    Molecule benzene_saved = fixture.benzene.get_source().lock();
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved, benzene_saved));
    }
    
    MolecularGraph loaded;
    Molecule benzene_loaded;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, benzene_loaded));
    }
    fixture.blank.imp = loaded;
    fixture.mol_fixture.blankmol.imp = benzene_loaded;
    
    check_eq(saved->NumVertices(), loaded->NumVertices());
    check_eq(saved->NumEdges(), loaded->NumEdges());
    check_eq(fixture.mol_fixture.blankmol.get_atms().size(),
             fixture.mol_fixture.benzene.get_atms().size());
    check_eq(fixture.mol_fixture.blankmol.get_bnds().size(),
             fixture.mol_fixture.benzene.get_bnds().size());
    
    for (Atom atm : fixture.mol_fixture.blankmol.get_atms())
      check(fixture.blank.get_g().HasVertex(loaded->GetVertex(atm).get()));
    for (Bond bnd : fixture.mol_fixture.blankmol.get_bnds())
      check(fixture.blank.get_g().HasEdge(loaded->GetEdge(bnd).get()));
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixmolg_serial, ixserial<IXMolecularGraph>);
*/
  
// ============================================================================
// == CONSTRUCTION ============================================================
// ============================================================================
  
  MGVertex::MGVertex() : _dat(nullptr) { }
  MGVertex::MGVertex(const MGVertex& v) : _dat(v._dat) { }
  MGVertex::MGVertex(MGVertex&& v) noexcept : _dat(std::move(v._dat)) { }
  MGVertex& MGVertex::operator=(const MGVertex &v) {
    if (&v != this) _dat = v._dat;
    return *this;
  }
  MGVertex& MGVertex::operator=(MGVertex &&v) {
    _dat = std::move(v._dat);
    return *this;
  }
  MGVertex::MGVertex(Atom& a, MolecularGraph& graph)
  : _dat(std::make_shared<MGVertexData>(a.weak_from_this(), graph.weak_from_this())) { }
  
  MGEdge::MGEdge() : _dat(nullptr) { }
  MGEdge::MGEdge(const MGEdge& e) : _dat(e._dat) { }
  MGEdge::MGEdge(MGEdge&& e) noexcept : _dat(std::move(e._dat)) { }
  MGEdge& MGEdge::operator=(const MGEdge& e) {
    if (&e != this) _dat = e._dat;
    return *this;
  }
  MGEdge& MGEdge::operator=(MGEdge &&e) {
    _dat = std::move(e._dat);
    return *this;
  }
  MGEdge::MGEdge(Bond& b, MolecularGraph& graph)
  : _dat(std::make_shared<MGEdgeData>(b.weak_from_this(), graph.weak_from_this())) { }
  
  MolecularGraph::MolecularGraph(Molecule& mol)
  : BaseGraph<MGVertex, MGEdge, sMolecularGraph>(), _at2v(), _bn2e(),
  _mol(mol.shared_from_this()), _subg() { }
  
/*  test_case("IXMolecularGraph construction") {
    using G = test::TestMolecularGraph;
    using V = test::TestMolecularVertex;
    using E = test::TestMolecularEdge;
    check_nothrow(G t(CreateMolecule()));
    check_nothrow(V t(test::CreateGenericTestAtom().imp,
                      test::CreateGenericTestMolecularGraph().imp));
    check_nothrow(E t(test::CreateGenericTestBond().imp,
                      test::CreateGenericTestMolecularGraph().imp));
  }
  
  test_case_fixture(test::MolecularGraphTestFixture, "IXMolecularGraph internals") {
    check(blank.get_at2v().empty());
    check(blank.get_bn2e().empty());
    check(blank.get_v().empty());
    check(blank.get_e().empty());
    check(blank.get_n().empty());
    check_eq(blank.get_source().lock(), mol_fixture.blankmol.imp);
  }
  */
  
// ============================================================================
// == FINDING ITEMS ===========================================================
// ============================================================================
  
/*  test_case_fixture(test::MolecularGraphTestFixture, "IXMolecularGraph has items") {
    for (Atom atm : mol_fixture.a_randmol) {
      check_false(benzene.HasVertex(atm));
      check(random.HasVertex(atm));
    }
    
    for (MGVertex v : random.get_v()) {
      check_false(benzene.HasVertex(v));
      check(random.HasVertex(v));
    }
    
    for (Bond bnd : mol_fixture.b_randmol) {
      check_false(benzene.HasEdge(bnd));
      check(random.HasEdge(bnd));
    }
    
    for (MGEdge e : random.get_e()) {
      check_false(benzene.HasEdge(e));
      check(random.HasEdge(e));
    }
    
    auto vertices = random.get_v();
    for (auto st : mol_fixture.e_randmol) {
      check_false(benzene.HasEdge(vertices[st.first], vertices[st.second]));
      check(random.HasEdge(vertices[st.first], vertices[st.second]));
      check(random.HasEdge(vertices[st.second], vertices[st.first]));
    }
    
    // Null checks
    check_false(benzene.HasVertex(Atom()));
    check_false(benzene.HasVertex(MGVertex()));
    check_false(benzene.HasEdge(Bond()));
    check_false(benzene.HasEdge(MGEdge()));
    check_false(benzene.HasEdge(MGVertex(), MGVertex()));
  }
 */
  
  
// ============================================================================
// == GETTING =================================================================
// ============================================================================
  
  Atom& MGVertex::GetAtom() const { return *_dat->atom.lock(); }
  MolecularGraph& MGVertex::GetGraph() const { return *_dat->graph.lock(); }
  
  Bond& MGEdge::GetBond() const { return *_dat->bond.lock(); }
  MolecularGraph& MGEdge::GetGraph() const { return *_dat->graph.lock(); }
  
  MGEdge MolecularGraph::GetEdge(Bond& bnd) const {
    if (!HasEdge(bnd)) throw std::out_of_range("No such edge");
    return _bn2e.at(bnd.shared_from_this());
  }
  
  MGVertex MolecularGraph::GetVertex(Atom& atm) const {
    if (!HasVertex(atm)) throw std::out_of_range("No such vertex");
    return _at2v.at(atm.shared_from_this());
  }
  
  Molecule& MolecularGraph::GetMolecule() const {
    if (IsSubgraph())
      throw std::runtime_error("Subgraphs do not relate to a molecule");
    return *_mol.lock();
  }
  
  MolecularGraph& MolecularGraph::GetSuperGraph() const {
    if (!IsSubgraph())
      throw std::runtime_error("Cannot get supergraph as not a subgraph");
    return *_subg.lock();
  }
  
/*  test_case_fixture(test::MolecularGraphTestFixture, "IXMolecularGraph getting") {
    // Numbers of things
    check_eq(mol_fixture.benzene.NumAtoms(), benzene.NumVertices());
    check_eq(mol_fixture.randmol.NumBonds(), random.NumEdges());
    
    // Degree of vertices
    auto vertices = random.get_v();
    for (size_t i = 0; i < random.NumVertices(); ++i) {
      check_eq(mol_fixture.a_randmol[i]->NumBonds(), random.Degree(vertices[i]));
      check_eq(std::numeric_limits<size_>::max(), benzene.Degree(vertices[i]));
    }
    
    // Getting vertices
    for (Atom atm : mol_fixture.a_randmol) {
      check_eq(vertices[atm->GetIndex()], random.GetVertex(atm));
      check_eq(MGVertex(), benzene.GetVertex(atm));
    }
    test::TestMolecularGraph::Vertices got_vertices(random.GetVertices().first,
                                                    random.GetVertices().second);
    check_eq(vertices, got_vertices);
    
    // Getting edges
    auto edges = random.get_e();
    for (Bond bnd : mol_fixture.b_randmol) {
      check_eq(edges[bnd->GetIndex()], random.GetEdge(bnd));
      check_eq(MGEdge(), benzene.GetEdge(bnd));
    }
    test::TestMolecularGraph::Edges got_Edges(random.GetEdges().first,
                                              random.GetEdges().second);
    check_eq(edges, got_Edges);
    
    for (size_t i = 0; i < mol_fixture.e_randmol.size(); ++i) {
      size_t u, v;
      std::tie(u, v) = mol_fixture.e_randmol[i];
      check_eq(edges[i], random.GetEdge(vertices[u], vertices[v]));
      check_eq(MGEdge(), benzene.GetEdge(vertices[u], vertices[v]));
      check_eq(vertices[u], random.GetSource(edges[i]));
      check_eq(vertices[v], random.GetTarget(edges[i]));
      check_eq(MGVertex(), benzene.GetSource(edges[i]));
      check_eq(MGVertex(), benzene.GetTarget(edges[i]));
      check_eq(std::make_pair(vertices[u], vertices[v]),
               random.GetVertices(edges[i]));
      check_eq(std::make_pair(MGVertex(), MGVertex()),
               benzene.GetVertices(edges[i]));
    }
    
    // checking neighbours
    std::set<size_t> expected, obtained;
    for (size_t i = 0; i < 12; ++i) {
      obtained.clear(); expected.clear();
      MGVertex v = benzene.GetVertex(mol_fixture.a_benzene[i]);
      auto nbrs = benzene.GetNeighbours(v);
      for (; nbrs.first != nbrs.second; ++nbrs.first)
        obtained.insert((*nbrs.first)->GetAtom()->GetIndex());
      for (auto ab : mol_fixture.e_benzene) {
        if (ab.first == i) expected.insert(ab.second);
        else if (ab.second == i) expected.insert(ab.first);
      }
      check_eq(expected, obtained);
    }
  }
 */
  
// ============================================================================
// == Subgraph Generation =====================================================
// ============================================================================
  sMolecularGraph MolecularGraph::Subgraph(std::vector<MGVertex> &verts) {
    sMolecularGraph G = std::make_shared<MolecularGraph>();
    G->_subg = weak_from_this();
    
    for (const MGVertex& v : verts) {
      if (!HasVertex(v)) throw std::runtime_error("Non-member vertex");
      G->_at2v.emplace(v.GetAtom().shared_from_this(), v);
      G->graph_type::AddVertex(v);
    }
    
    for (const MGEdge& e : GetEdges()) {
      MGVertex u = GetSourceVertex(e);
      MGVertex v = GetTargetVertex(e);
      if (!G->HasVertex(u) || !G->HasVertex(v)) continue;
      G->_bn2e.emplace(e.GetBond().shared_from_this(), e);
      G->graph_type::AddEdge(u, v, e);
    }
    return G;
  }
  
  sMolecularGraph MolecularGraph::Subgraph(std::vector<MGVertex> &verts,
                                           std::vector<MGEdge> &edges) {
    sMolecularGraph G = std::make_shared<MolecularGraph>();
    G->_subg = weak_from_this();
    
    for (const MGVertex& v : verts) {
      if (!HasVertex(v)) throw std::runtime_error("Non-member vertex");
      G->_at2v.emplace(v.GetAtom().shared_from_this(), v);
      G->graph_type::AddVertex(v);
    }
    
    for (const MGEdge& e : edges) {
      if (!HasEdge(e)) throw std::runtime_error("Non-member edge");
      MGVertex u = GetSourceVertex(e);
      MGVertex v = GetTargetVertex(e);
      if (!G->HasVertex(u) || !G->HasVertex(v)) continue;
      G->_bn2e.emplace(e.GetBond().shared_from_this(), e);
      G->graph_type::AddEdge(u, v, e);
    }
    return G;
  }
  
  
// ============================================================================
// == STRUCTURE MODIFICATION ==================================================
// ============================================================================
  
  void MolecularGraph::Clear() {
    _at2v.clear();
    _bn2e.clear();
    graph_type::Clear();
    _mol.reset();
    _cond->Clear();
    _cond.reset();
  }
  
  MGVertex MolecularGraph::AddVertex(Atom& atm) {
    MGVertex v(atm, *this);
    _at2v.emplace(atm.shared_from_this(), v);
    graph_type::AddVertex(v);
    return v;
  }
  
  void MolecularGraph::RemoveVertex(MGVertex& v) {
    _at2v.erase(v.GetAtom().shared_from_this());
    for (const MGVertex& u : graph_type::GetNeighbours(v)) {
      MGEdge e = GetEdge(u, v);
      _bn2e.erase(e.GetBond().shared_from_this());
    }
    graph_type::RemoveVertex(v);
  }
  
 /* test_case_fixture(test::MolecularGraphTestFixture,
                    "IXMolecularGraph adding and removing vertices") {
    
    // New Vertex
    Atom atm1 = test::CreateGenericTestAtom().imp;
    MGVertex v1 = blank.AddVertex(atm1);
    check_ne(MGVertex(), v1);
    check_eq(1, blank.NumVertices());
    check_eq(1, blank.get_v().size());
    check_eq(1, blank.get_n().size());
    check_eq(1, blank.get_at2v().size());
    check_eq(0, blank.get_n().at(v1).size());
    
    // Remove vertex
    check_nothrow(blank.RemoveVertex(v1));
    check_eq(0, blank.NumVertices());
    check_eq(0, blank.get_v().size());
    check_eq(0, blank.get_n().size());
    check_eq(0, blank.get_at2v().size());
    check_false(blank.HasVertex(atm1));
    check_false(blank.HasVertex(v1));
    
    // No null cases as normal usage is via IXMolecule
    
    // Removing a vertex should remove all edges around it as well
    check_nothrow(benzene.RemoveVertex(benzene.get_v().front()));
    check_eq(11, benzene.NumVertices());
    check_eq(9, benzene.NumEdges());
    check_eq(11, benzene.get_v().size());
    check_eq(11, benzene.get_n().size());
    check_eq(11, benzene.get_at2v().size());
    check_eq(2, benzene.get_n().at(benzene.get_v()[0]).size());
    check_eq(2, benzene.get_n().at(benzene.get_v()[4]).size());
    check_eq(0, benzene.get_n().at(benzene.get_v()[5]).size());
  }
  */
  
  MGEdge MolecularGraph::AddEdge(Bond& bnd) {
    Atom& source = bnd.GetSourceAtom();
    Atom& target = bnd.GetTargetAtom();
    // Edges will only be added when a bond is added to the molecule
    MGVertex u = GetVertex(source);
    MGVertex v = GetVertex(target);
    MGEdge e(bnd, *this);
    _bn2e.emplace(bnd.shared_from_this(), e);
    graph_type::AddEdge(u, v, e);
    return e;
  }
  
  void MolecularGraph::RemoveEdge(MGEdge& e) {
    _bn2e.erase(e.GetBond().shared_from_this());
    graph_type::RemoveEdge(e);
  }
  
  void MolecularGraph::RemoveEdge(MGVertex& u, MGVertex& v) {
    MGEdge e = GetEdge(u, v);
    RemoveEdge(e);
  }
  
  bool MolecularGraph::HasVertex(Atom &v) const {
    return _at2v.find(v.shared_from_this()) != _at2v.end();
  }
  
  bool MolecularGraph::HasEdge(Bond &e) const {
    return _bn2e.find(e.shared_from_this()) != _bn2e.end(); }
  
/*  test_case_fixture(test::MolecularGraphTestFixture,
                    "IXMolecularGraph adding and removing edges") {
    Atom a1 = test::CreateGenericTestAtom().imp;
    Atom a2 = test::CreateGenericTestAtom().imp;
    Atom a3 = test::CreateGenericTestAtom().imp;
    Atom a4 = test::CreateGenericTestAtom().imp;
    Bond b1 = test::CreateGenericTestBond(a1,a2).imp;
    Bond b2 = test::CreateGenericTestBond(a4,a3).imp;
    MGVertex v1 = blank.AddVertex(a1);
    MGVertex v2 = blank.AddVertex(a2);
    
    // New edge
    MGEdge e1 = blank.AddEdge(b1);
    check_ne(MGEdge(), e1);
    check_eq(1, blank.NumEdges());
    check_eq(1, blank.get_e().size());
    check_eq(2, blank.get_n().size());
    check_eq(1, blank.get_bn2e().size());
    check_eq(1, blank.get_n().at(v1).size());
    check_eq(1, blank.get_n().at(v2).size());
    
    // New edge without vertices should add vertices
    MGEdge e2 = blank.AddEdge(b2);
    check_ne(MGEdge(), e2);
    check(blank.HasVertex(a3));
    check(blank.HasVertex(a4));
    check_eq(2, blank.NumEdges());
    check_eq(2, blank.get_e().size());
    check_eq(4, blank.get_n().size());
    check_eq(2, blank.get_bn2e().size());
    
    // Remove an edge
    check_nothrow(benzene.RemoveEdge(benzene.get_e().front()));
    check_eq(12, benzene.NumVertices());
    check_eq(11, benzene.NumEdges());
    check_eq(12, benzene.get_v().size());
    check_eq(12, benzene.get_n().size());
    check_eq(12, benzene.get_at2v().size());
    check_eq(11, benzene.get_bn2e().size());
    check_eq(2, benzene.get_n().at(benzene.get_v()[0]).size());
    check_eq(2, benzene.get_n().at(benzene.get_v()[1]).size());
    
    // Remove an edge between vertices
    check_nothrow(blank.RemoveEdge(v2, v1));
    check_eq(1, blank.NumEdges());
    check_eq(1, blank.get_e().size());
    check_eq(4, blank.get_n().size());
    check_eq(1, blank.get_bn2e().size());
  }
*/
  
//  std::pair<IXMolecularGraph::CompIter, IXMolecularGraph::CompIter>
//  IXMolecularGraph::GetConnectedComponents() {
//    _c.clear();
//    std::vector<std::vector<IXMGVertex*>> ptr_components;
//    size_ num = _g.ConnectedComponents(ptr_components);
//    _c.reserve(num);
//    for (auto component : ptr_components) {
//      std::vector<MGVertex> c;
//      c.reserve(component.size());
//      for (IXMGVertex* v : component)
//        c.emplace_back(v->shared_from_this());
//      _c.emplace_back(c.begin(), c.end());
//    }
//    return {_c.cbegin(), _c.cend()};
//  }
  
  CondensedMolecularGraph& MolecularGraph::GetCondensedGraph() {
    if (!IsFrozen())
      throw std::runtime_error("Can only condense a frozen graph");
    if (!_cond) _cond = Condense(*this);
    return *_cond;
  }
  
  test_suite_close();
}


