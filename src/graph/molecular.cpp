#include <algorithm>
#include <limits>

#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/serialise.hpp>

#include <indigox/utils/doctest_proxy.hpp>
#include <indigox/test/molecular_test.hpp>
#include <indigox/test/atom_test.hpp>
#include <indigox/test/bond_test.hpp>

namespace indigox::graph {
  test_suite_open("IXMolecularGraph");
  
// ============================================================================
// == SERIALISATION ===========================================================
// ============================================================================
  
  template <typename Archive>
  void IXMGVertex::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("source_atom", _atom));
  }
  
  template <typename Archive>
  void IXMGVertex::load_and_construct(Archive &archive,
                                      cereal::construct<IXMGVertex> &construct,
                                      const uint32_t) {
    Atom atom;
    archive(INDIGOX_SERIAL_NVP("source_atom", atom));
    construct(atom);
  }
  INDIGOX_SERIALISE_SPLIT(IXMGVertex);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMGVertex serialisation", T, ixmolv_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestMolecularVertex saved(test::CreateGenericTestAtom().imp);
  
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp, saved.GetAtom()));
    }
    
    Atom atom_loaded;
    graph::MGVertex loaded;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, atom_loaded));
    }
    saved.imp = loaded;
    
    check_eq(loaded->GetAtom(), atom_loaded);
    check_eq(saved.get_atom().lock(), atom_loaded);
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixmolv_serial, ixserial<IXMGVertex>);
  
  template <typename Archive>
  void IXMGEdge::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("source_bond", _bond));
  }
  
  template <typename Archive>
  void IXMGEdge::load_and_construct(Archive &archive,
                                    cereal::construct<IXMGEdge> &construct,
                                    const uint32_t) {
    Bond bond;
    archive(INDIGOX_SERIAL_NVP("source_bond", bond));
    construct(bond);
  }
  INDIGOX_SERIALISE_SPLIT(IXMGEdge);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMGEdge serialisation", T, ixmole_serial) {
    using In = typename T::t1;
    using Out = typename cereal::traits::detail::get_output_from_input<In>::type;
    
    test::TestMolecularEdge saved(test::CreateGenericTestBond().imp);
    
    std::ostringstream os;
    {
      Out oar(os);
      check_nothrow(oar(saved.imp, saved.GetBond()));
    }
    
    Bond bond_loaded;
    graph::MGEdge loaded;
    std::istringstream is(os.str());
    {
      In iar(is);
      check_nothrow(iar(loaded, bond_loaded));
    }
    saved.imp = loaded;
    
    check_eq(loaded->GetBond(), bond_loaded);
    check_eq(saved.get_bond().lock(), bond_loaded);
  }
  DOCTEST_TEST_CASE_TEMPLATE_INSTANTIATE(ixmole_serial, ixserial<IXMGEdge>);
  
  template <typename Archive>
  void IXMolecularGraph::save(Archive &archive, const uint32_t) const {
    archive(INDIGOX_SERIAL_NVP("source_molecule", _source),
            INDIGOX_SERIAL_NVP("vertices", _v),
            INDIGOX_SERIAL_NVP("edges", _e));
  }
  
  template <typename Archive>
  void IXMolecularGraph::load_and_construct(Archive &archive,
                                            cereal::construct<IXMolecularGraph> &construct,
                                            const uint32_t) {
    // Reload everything
    Molecule mol;
    archive(INDIGOX_SERIAL_NVP("source_molecule", mol));
    construct(mol);
    archive(INDIGOX_SERIAL_NVP("vertices", construct->_v),
            INDIGOX_SERIAL_NVP("edges", construct->_e));
    
    // Rebuild the underlying graph
    for (MGVertex v : construct->_v) {
      construct->_at2v.emplace(v->GetAtom(), v);
      construct->_n.emplace(v, VertContain());
      construct->_g.AddVertex(v.get());
    }
    for (MGEdge e : construct->_e) {
      Bond b = e->GetBond();
      construct->_bn2e.emplace(b, e);
      MGVertex u = construct->_at2v.at(b->GetSourceAtom());
      MGVertex v = construct->_at2v.at(b->GetTargetAtom());
      construct->_n[u].emplace_back(v);
      construct->_n[v].emplace_back(u);
      construct->_g.AddEdge(u.get(), v.get(), e.get());
    }
  }
  INDIGOX_SERIALISE_SPLIT(IXMolecularGraph);
  
  DOCTEST_TEST_CASE_TEMPLATE_DEFINE("IXMolecularGraph serialisation", T, ixmolg_serial) {
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
  
// ============================================================================
// == CONSTRUCTION ============================================================
// ============================================================================
  
  test_case_fixture(test::MolecularGraphTestFixture, "IXMolecularGraph construction") {
    check_nothrow(test::TestMolecularGraph test1(CreateMolecule()));
    check_nothrow(test::TestMolecularVertex(test::CreateGenericTestAtom().imp));
    check_nothrow(test::TestMolecularEdge(test::CreateGenericTestBond().imp));
    
    check(blank.get_at2v().empty());
    check(blank.get_bn2e().empty());
    check(blank.get_v().empty());
    check(blank.get_e().empty());
    check(blank.get_n().empty());
    check_eq(blank.get_source().lock(), mol_fixture.blankmol.imp);
  }
  
// ============================================================================
// == FINDING ITEMS ===========================================================
// ============================================================================
  
  test_case_fixture(test::MolecularGraphTestFixture, "IXMolecularGraph has items") {
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
  
  
// ============================================================================
// == GETTING =================================================================
// ============================================================================
  
  size_ IXMolecularGraph::Degree(const MGVertex v) const {
    if (!HasVertex(v)) return std::numeric_limits<size_>::max();
    return _g.Degree(v.get());
  }
  
  MGEdge IXMolecularGraph::GetEdge(const MGVertex u, const MGVertex v) const {
    if (!HasEdge(u, v)) return MGEdge();
    return _g.GetEdge(u.get(), v.get())->shared_from_this();
  }
  
  MGEdge IXMolecularGraph::GetEdge(const Bond bnd) const {
    if (!HasEdge(bnd)) return MGEdge();
    return _bn2e.at(bnd);
  }
  
  MGVertex IXMolecularGraph::GetVertex(const Atom atm) const {
    if (!HasVertex(atm)) return MGVertex();
    return _at2v.at(atm);
  }
  
  MGVertex IXMolecularGraph::GetSource(const MGEdge e) const {
    if (!HasEdge(e)) return MGVertex();
    return _g.GetSource(e.get())->shared_from_this();
  }
  
  MGVertex IXMolecularGraph::GetTarget(const MGEdge e) const {
    if (!HasEdge(e)) return MGVertex();
    return _g.GetTarget(e.get())->shared_from_this();
  }
  
  std::pair<MGVertex, MGVertex>
  IXMolecularGraph::GetVertices(const MGEdge e) const {
    if (!HasEdge(e)) return {MGVertex(), MGVertex()};
    auto res = _g.GetVertices(e.get());
    return {res.first->shared_from_this(), res.second->shared_from_this()};
  }
  
  test_case_fixture(test::MolecularGraphTestFixture, "IXMolecularGraph getting") {
    // Numbers of things
    check_eq(mol_fixture.benzene.NumAtoms(), benzene.NumVertices());
    check_eq(mol_fixture.randmol.NumBonds(), random.NumEdges());
    
    // Degree of vertices
    auto vertices = random.get_v();
    for (size_ i = 0; i < random.NumVertices(); ++i) {
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
    
    for (size_ i = 0; i < mol_fixture.e_randmol.size(); ++i) {
      size_ u, v;
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
    std::set<size_> expected, obtained;
    for (size_ i = 0; i < 12; ++i) {
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
  
// ============================================================================
// == STRUCTURE MODIFICATION ==================================================
// ============================================================================
  
  MGVertex IXMolecularGraph::AddVertex(const Atom atm) {
    MGVertex v = MGVertex(new IXMGVertex(atm));
    _at2v.emplace(atm, v);
    _v.emplace_back(v);
    _n.emplace(v, VertContain());
    _g.AddVertex(v.get());
    return v;
  }
  
  void IXMolecularGraph::RemoveVertex(const MGVertex v) {
    std::vector<MGVertex> nbrs(_n[v].begin(), _n[v].end());
    for (MGVertex n : nbrs) RemoveEdge(v, n);
    _g.RemoveVertex(v.get());
    _at2v.erase(v->GetAtom());
    _v.erase(std::find(_v.begin(), _v.end(), v));
    _n.erase(v);
  }
  
  test_case_fixture(test::MolecularGraphTestFixture,
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
  
  MGEdge IXMolecularGraph::AddEdge(const Bond bnd) {
    Atom source = bnd->GetSourceAtom();
    Atom target = bnd->GetTargetAtom();
    if (!HasVertex(source)) AddVertex(source);
    if (!HasVertex(target)) AddVertex(target);
    // Edges will only be added when a bond is added to the molecule
    MGVertex u = GetVertex(source);
    MGVertex v = GetVertex(target);
    MGEdge e = MGEdge(new IXMGEdge(bnd));
    _bn2e.emplace(bnd, e);
    _e.emplace_back(e);
    _n[u].emplace_back(v);
    _n[v].emplace_back(u);
    _g.AddEdge(u.get(), v.get(), e.get());
    return e;
  }
  
  void IXMolecularGraph::RemoveEdge(const MGEdge e) {
    _g.RemoveEdge(e.get());
    _bn2e.erase(e->GetBond());
    _e.erase(std::find(_e.begin(), _e.end(), e));
    MGVertex u = _at2v[e->GetBond()->GetSourceAtom()];
    MGVertex v = _at2v[e->GetBond()->GetTargetAtom()];
    _n[u].erase(std::find(_n[u].begin(), _n[u].end(), v));
    _n[v].erase(std::find(_n[v].begin(), _n[v].end(), u));
  }
  
  void IXMolecularGraph::RemoveEdge(const MGVertex u, const MGVertex v) {
    RemoveEdge(GetEdge(u, v));
  }
  
  test_case_fixture(test::MolecularGraphTestFixture,
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
  
  void IXMolecularGraph::Clear() {
    _source.reset();
    _g.Clear();
    _at2v.clear();
    _bn2e.clear();
    _v.clear();
    _e.clear();
    _n.clear();
  }
  
  test_case_fixture(test::MolecularGraphTestFixture,
                    "IXMolecularGraph clearing methods") {
    check_false(benzene.get_source().expired());
    check_ne(0, benzene.get_g().NumVertices());
    check_false(benzene.get_at2v().empty());
    check_false(benzene.get_bn2e().empty());
    check_false(benzene.get_v().empty());
    check_false(benzene.get_e().empty());
    check_false(benzene.get_n().empty());
    
    check_nothrow(benzene.Clear());
    
    check(benzene.get_source().expired());
    check_eq(0, benzene.get_g().NumVertices());
    check(benzene.get_at2v().empty());
    check(benzene.get_bn2e().empty());
    check(benzene.get_v().empty());
    check(benzene.get_e().empty());
    check(benzene.get_n().empty());
  }
  
  
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
  
  test_suite_close();
}

