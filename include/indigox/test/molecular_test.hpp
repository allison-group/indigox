#ifndef INDIGOX_TEST_MOLECULAR_TEST_HPP
#define INDIGOX_TEST_MOLECULAR_TEST_HPP

#include "molecule_test.hpp"

namespace indigox::test {
  struct TestMolecularVertex {
//    graph::MGVertex imp;
//
//    // Private wrapping functions
//    TestMolecularVertex() = delete;
//    TestMolecularVertex(Atom a, graph::MolecularGraph g)
//    : imp(new graph::IXMGVertex(a, g)) { }
//
//    // Public wrapping functions
//    Atom GetAtom() const { return imp->GetAtom(); }
//    graph::MolecularGraph GetGraph() const { return imp->GetGraph(); }
//
//    // Internals access
//    _Atom get_atom() const { return imp->_atom; }
//    graph::_MolecularGraph get_graph() const { return imp->_graph; }
//  };
//
//  struct TestMolecularEdge {
//    graph::MGEdge imp;
//
//    // Private wrapping functions
//    TestMolecularEdge() = delete;
//    TestMolecularEdge(Bond b, graph::MolecularGraph g)
//    : imp(new graph::IXMGEdge(b, g)) { }
//
//    // Public wrapping functions
//    Bond GetBond() const { return imp->GetBond(); }
//    graph::MolecularGraph GetGraph() const { return imp->GetGraph(); }
//
//    // Internals access
//    _Bond get_bond() const { return imp->_bond; }
//    graph::_MolecularGraph get_graph() const { return imp->_graph; }
  };
  
  struct TestMolecularGraph {
    // Typedefs
//    using Vertices = indigox::graph::IXMolecularGraph::VertContain;
//    using Edges = indigox::graph::IXMolecularGraph::EdgeContain;
//    using AllNeighbours = indigox::graph::IXMolecularGraph::NbrsContain;
//    using MGEdge = indigox::graph::MGEdge;
//    using MGVertex = indigox::graph::MGVertex;
//    using EdgeIter = indigox::graph::IXMolecularGraph::EdgeIter;
//    using NbrsIter = indigox::graph::IXMolecularGraph::NbrsIter;
//    using VertIter = indigox::graph::IXMolecularGraph::VertIter;
//    using GraphType = indigox::graph::IXMolecularGraph::graph_type;
//    using AtomMap = indigox::graph::IXMolecularGraph::AtomMap;
//    using BondMap = indigox::graph::IXMolecularGraph::BondMap;
//
//    indigox::graph::MolecularGraph imp;
//
//    // Private wrapping functions
//    TestMolecularGraph() = default;
//    TestMolecularGraph(const Molecule& mol) : imp(new graph::IXMolecularGraph(mol)) { }
//    MGEdge AddEdge(const Bond& b) { return imp->AddEdge(b); }
//    MGVertex AddVertex(const Atom& a) { return imp->AddVertex(a); }
//    void Clear() { imp->Clear(); }
//    void RemoveEdge(const MGEdge& e) { imp->RemoveEdge(e); }
//    void RemoveEdge(const MGVertex& u, const MGVertex& v) { imp->RemoveEdge(u,v); }
//    void RemoveVertex(const MGVertex& v) { imp->RemoveVertex(v); }
//
//    // Public wrapping functions
//    size_t Degree(const MGVertex& v) const { return imp->Degree(v); }
//    MGEdge GetEdge(const MGVertex& u, const MGVertex& v) const { return imp->GetEdge(u,v); }
//    MGEdge GetEdge(const Bond& b) const { return imp->GetEdge(b); }
//    MGVertex GetVertex(const Atom& a) const { return imp->GetVertex(a); }
//    std::pair<EdgeIter,EdgeIter> GetEdges() const { return imp->GetEdges(); }
//    std::pair<NbrsIter,NbrsIter> GetNeighbours(const MGVertex& v) const { return imp->GetNeighbours(v); }
//    MGVertex GetSource(const MGEdge& e) const { return imp->GetSource(e); }
//    MGVertex GetTarget(const MGEdge& e) const { return imp->GetTarget(e); }
//    std::pair<MGVertex,MGVertex> GetVertices(const MGEdge& e) const { return imp->GetVertices(e); }
//    std::pair<VertIter,VertIter> GetVertices() const { return imp->GetVertices(); }
//    bool HasVertex(const Atom& v) const { return imp->HasVertex(v); }
//    bool HasVertex(const MGVertex& v) const { return imp->HasVertex(v); }
//    bool HasEdge(const Bond& e) const { return imp->HasEdge(e); }
//    bool HasEdge(const MGEdge& e) const { return imp->HasEdge(e); }
//    bool HasEdge(const MGVertex& u, const MGVertex& v) const { return imp->HasEdge(u,v); }
//    size_t NumEdges() const { return imp->NumEdges(); }
//    size_t NumVertices() const { return imp->NumVertices(); }
//
//    // Internals access
//    _Molecule get_source() const { return imp->_source; }
//    const GraphType& get_g() const { return imp->_g; }
//    const AtomMap& get_at2v() const { return imp->_at2v; }
//    const BondMap& get_bn2e() const { return imp->_bn2e; }
//    const Vertices& get_v() const { return imp->_v; }
//    const Edges& get_e() const { return imp->_e; }
//    const AllNeighbours& get_n() const { return imp->_n; }
  };

//  inline TestMolecularGraph CreateGenericTestMolecularGraph() {
//    return TestMolecularGraph(CreateMolecule());
//  }
  
  struct MolecularGraphTestFixture {
//    TestMolecularGraph benzene, random, blank;
//    MoleculeTestFixture mol_fixture;
//    
//    MolecularGraphTestFixture() {
//      mol_fixture.BuildBenzene();
//      mol_fixture.BuildRandomMolecule();
//      mol_fixture.blankmol.Init();
//      benzene.imp = mol_fixture.benzene.get_g();
//      random.imp = mol_fixture.randmol.get_g();
//      blank.imp = mol_fixture.blankmol.get_g();
//    }
  };
}

#endif /* INDIGOX_TEST_MOLECULAR_TEST_HPP */
