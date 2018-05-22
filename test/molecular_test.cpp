#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include "class_test_wrappers.hpp"

namespace indigox::test {
  struct MolecularGraphFixture {
    Molecule mol;
    indigox::test::IXMolecularGraph G;
    size_ num_atoms, num_bonds;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<graph::MGVertex> verts;
    std::vector<graph::MGEdge> edges;
    std::vector<std::pair<size_, size_>> incident_edges;
    std::pair<size_, size_> incident_fail_edge;
    std::map<size_, size_> expected_degrees;
    std::map<graph::MGVertex, size_> vert_ids;
    MolecularGraphFixture() : mol(CreateMolecule()), G(mol) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<size_> at_dis(15,35);
      std::uniform_int_distribution<size_> bn_dis(35,55);
      
      // Random number of atoms
      num_atoms = at_dis(gen);
      atoms.reserve(num_atoms); verts.reserve(num_atoms);
      for (int i = 0; i < num_atoms; ++i) {
        expected_degrees.emplace(i,0);
        atoms.emplace_back(indigox::test::IXAtom::GetNewAtom(mol));
        verts.emplace_back(G.AddVertex(atoms.back()));
        vert_ids.emplace(verts.back(), i);
      }
      
      // Random number of bonds
      num_bonds = bn_dis(gen);
      std::vector<std::pair<size_, size_>> possible_edges;
      possible_edges.reserve(num_atoms * (num_atoms + 1) / 2);
      bonds.reserve(num_bonds); edges.reserve(num_bonds); incident_edges.reserve(num_bonds);
      for (size_ i = 0; i < num_atoms; ++i) {
        for (size_ j = i + 1; j < num_atoms; ++j ) possible_edges.emplace_back(i,j);
      }
      std::sample(possible_edges.begin(), possible_edges.end(),
                  std::back_inserter(incident_edges), num_bonds + 1, gen);
      incident_fail_edge = incident_edges.back();
      incident_edges.resize(num_bonds);
      for (auto& i : incident_edges) {
        expected_degrees.at(i.first)++;
        expected_degrees.at(i.second)++;
        bonds.emplace_back(indigox::test::IXBond::GetNewBond(atoms[i.first], atoms[i.second], mol));
        edges.emplace_back(G.AddEdge(bonds.back()));
      }
    }
  };
}

using namespace indigox;
using namespace indigox::graph;

namespace std {
  std::ostream& operator<<(std::ostream& os, const std::pair<indigox::graph::MGVertex, indigox::graph::MGVertex> p) {
    return (os << &p.first << "," << &p.second);
  }
}

BOOST_FIXTURE_TEST_SUITE(ixmolecular_graph, indigox::test::MolecularGraphFixture);
//BOOST_AUTO_TEST_SUITE(ixmolecular_graph);

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK(G.NumEdges() == num_bonds);
  BOOST_CHECK(G.NumVertices() == num_atoms);
  auto its_verts = G.GetVertices();
  auto its_edges = G.GetEdges();
  BOOST_CHECK(its_verts.second == (its_verts.first + num_atoms));
  BOOST_CHECK(its_edges.second == (its_edges.first + num_bonds));
  
  // Check that edges and vertices are as expected
  auto vert_itr = G.GetVertices();
  auto edge_itr = G.GetEdges();
  std::set<MGVertex> expect_vert(verts.begin(), verts.end());
  std::set<MGVertex> obtain_vert(vert_itr.first, vert_itr.second);
  std::set<MGEdge> expect_edge(edges.begin(), edges.end());
  std::set<MGEdge> obtain_edge(edge_itr.first, edge_itr.second);
  BOOST_CHECK_EQUAL_COLLECTIONS(expect_vert.begin(), expect_vert.end(),
                                obtain_vert.begin(), obtain_vert.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(expect_edge.begin(), expect_edge.end(),
                                obtain_edge.begin(), obtain_edge.end());
}

BOOST_AUTO_TEST_CASE(vertex_add_remove) {
  // Check constructions
  std::set<MGVertex> expect_verts(verts.begin(), verts.end());
  auto verts_its = G.GetVertices();
  std::set<MGVertex> built_verts(verts_its.first, verts_its.second);
  BOOST_CHECK_EQUAL_COLLECTIONS(expect_verts.begin(), expect_verts.end(),
                                built_verts.begin(), built_verts.end());
  
  // Remove the first vertices 5 and 7, updating the expected_degrees
  size_ removed = 0;
  for (auto nbrs = G.GetNeighbours(verts[5]); nbrs.first != nbrs.second; ++nbrs.first) {
    size_ id = vert_ids.at(*nbrs.first);
    expected_degrees.at(id)--;
    ++removed;
  }
  expected_degrees.at(vert_ids.at(verts[5])) = std::numeric_limits<size_>::max();
  G.RemoveVertex(verts[5]);
  for (auto nbrs = G.GetNeighbours(verts[7]); nbrs.first != nbrs.second; ++nbrs.first) {
    size_ id = vert_ids.at(*nbrs.first);
    expected_degrees.at(id)--;
    ++removed;
  }
  expected_degrees.at(vert_ids.at(verts[7])) = std::numeric_limits<size_>::max();
  G.RemoveVertex(verts[7]);
  num_atoms -= 2;
  num_bonds -= removed;
  
  // Check degrees
  std::vector<size_> built_degree, expect_degree;
  for (auto v : verts) built_degree.push_back(G.Degree(v));
  for (auto ed : expected_degrees) expect_degree.push_back(ed.second);
  BOOST_CHECK_EQUAL_COLLECTIONS(expect_degree.begin(), expect_degree.end(),
                                built_degree.begin(), built_degree.end());
  BOOST_CHECK(G.NumEdges() == num_bonds);
  BOOST_CHECK(G.NumVertices() == num_atoms);
}

BOOST_AUTO_TEST_CASE(edge_add_remove) {
  // Check constructions
  std::set<MGEdge> expect_edges(edges.begin(), edges.end());
  auto edges_its = G.GetEdges();
  std::set<MGEdge> built_edges(edges_its.first, edges_its.second);
  BOOST_CHECK_EQUAL_COLLECTIONS(expect_edges.begin(), expect_edges.end(),
                                built_edges.begin(), built_edges.end());
  
  // Add an edge where the vertices haven't been added
  atoms.emplace_back(indigox::test::IXAtom::GetNewAtom(mol));
  atoms.emplace_back(indigox::test::IXAtom::GetNewAtom(mol));
  bonds.emplace_back(indigox::test::IXBond::GetNewBond(atoms[atoms.size() - 1], atoms[atoms.size() - 2], mol));
  edges.emplace_back(G.AddEdge(bonds.back()));
  BOOST_CHECK(G.NumVertices() == num_atoms + 2);
  BOOST_CHECK(G.NumEdges() == num_bonds + 1);
  
  // Remove the first edge
  G.RemoveEdge(edges[0]);
  BOOST_CHECK(G.NumVertices() == num_atoms + 2);
  BOOST_CHECK(G.NumEdges() == num_bonds);
  
  // Remove the 7th edge by vertices
  G.RemoveEdge(verts[incident_edges[7].first], verts[incident_edges[7].second]);
  BOOST_CHECK(G.NumVertices() == num_atoms + 2);
  BOOST_CHECK(G.NumEdges() == num_bonds - 1);
}

BOOST_AUTO_TEST_CASE(has_vertex_edge) {
  // Check that all verts  and edges correctly have
  std::vector<bool> tests_verts(num_atoms, true);
  std::vector<bool> tests_edges(num_bonds, true);
  std::vector<bool> fill_verts, fill_edges;
  for (auto v : verts) fill_verts.push_back(G.HasVertex(v));
  for (auto e : edges) fill_edges.push_back(G.HasEdge(e));
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_verts.begin(), tests_verts.end(),
                                fill_verts.begin(), fill_verts.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_edges.begin(), tests_edges.end(),
                                fill_edges.begin(), fill_edges.end());
  
  // Check that all atoms and bonds correctly have
  fill_verts.clear(); fill_edges.clear();
  for (auto v : atoms) fill_verts.push_back(G.HasVertex(v));
  for (auto e : bonds) fill_edges.push_back(G.HasEdge(e));
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_verts.begin(), tests_verts.end(),
                                fill_verts.begin(), fill_verts.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_edges.begin(), tests_edges.end(),
                                fill_edges.begin(), fill_edges.end());
  
  // Check that all incident_Edges corrcetly have by vertex
  fill_edges.clear();
  for (auto uv : incident_edges) fill_edges.push_back(G.HasEdge(verts[uv.first], verts[uv.second]));
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_edges.begin(), tests_edges.end(),
                                fill_edges.begin(), fill_edges.end());
  
  // Check that has methods correctly fail
  Atom fail_atom = Atom(indigox::test::IXAtom::GetNewAtom(mol));
  MGVertex fail_vert = G.AddVertex(fail_atom); G.RemoveVertex(fail_vert);
  BOOST_CHECK(!G.HasVertex(Atom()));     // fail on null pass
  BOOST_CHECK(!G.HasVertex(fail_atom));  // fail on not part pass
  BOOST_CHECK(!G.HasVertex(MGVertex())); // fail on null pass
  BOOST_CHECK(!G.HasVertex(fail_vert));  // fail on not part pass
  
  Bond fail_bond = indigox::test::IXBond::GetNewBond(atoms[incident_fail_edge.first], atoms[incident_fail_edge.second], mol);
  MGEdge fail_edge = G.AddEdge(fail_bond); G.RemoveEdge(fail_edge);
  BOOST_CHECK(!G.HasEdge(fail_bond));
  BOOST_CHECK(!G.HasEdge(fail_edge));
  BOOST_CHECK(!G.HasEdge(Bond()));
  BOOST_CHECK(!G.HasEdge(MGEdge()));
  BOOST_CHECK(!G.HasEdge(fail_vert, verts[0]));
  BOOST_CHECK(!G.HasEdge(verts[0], verts[0]));
}

BOOST_AUTO_TEST_CASE(edge_vertex_get) {
  // Check that all verts and edges correctly get by atom/bond
  std::vector<MGVertex> tests_verts;
  std::vector<MGEdge> tests_edges;
  for (auto v : atoms) tests_verts.push_back(G.GetVertex(v));
  for (auto e : bonds) tests_edges.push_back(G.GetEdge(e));
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_verts.begin(), tests_verts.end(),
                                verts.begin(), verts.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_edges.begin(), tests_edges.end(),
                                edges.begin(), edges.end());
  
  // Check that all edges correctly get by vert pair
  tests_edges.clear();
  for (auto uv : incident_edges) tests_edges.push_back(G.GetEdge(verts[uv.first], verts[uv.second]));
  BOOST_CHECK_EQUAL_COLLECTIONS(tests_edges.begin(), tests_edges.end(),
                                edges.begin(), edges.end());
  
  
  // Check that get methods correctly fail (return null)
  Atom fail_atom = Atom(indigox::test::IXAtom::GetNewAtom(mol));
  MGVertex fail_vert = G.AddVertex(fail_atom); G.RemoveVertex(fail_vert);
  BOOST_CHECK(G.GetVertex(Atom()) == MGVertex());     // fail on null pass
  BOOST_CHECK(G.GetVertex(fail_atom) == MGVertex());  // fail on not part pass
  
  Bond fail_bond = indigox::test::IXBond::GetNewBond(atoms[incident_fail_edge.first], atoms[incident_fail_edge.second], mol);
  BOOST_CHECK(G.GetEdge(fail_bond) == MGEdge());
  BOOST_CHECK(G.GetEdge(Bond()) == MGEdge());
  BOOST_CHECK(G.GetEdge(fail_vert, verts[0]) == MGEdge());
  BOOST_CHECK(G.GetEdge(verts[0], verts[0]) == MGEdge());
  
  // Check that source/target get correctly
  std::vector<MGVertex> expected_st;
  tests_verts.clear();
  for (auto uv : incident_edges) {
    expected_st.push_back(verts[uv.first]);
    expected_st.push_back(verts[uv.second]);
  }
  for (auto e: edges) {
    tests_verts.push_back(G.GetSource(e));
    tests_verts.push_back(G.GetTarget(e));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_st.begin(), expected_st.end(),
                                tests_verts.begin(), tests_verts.end());
  
  MGEdge fail_edge = G.AddEdge(fail_bond); G.RemoveEdge(fail_edge);
  // Check that source/target fail correcly
  BOOST_CHECK(G.GetSource(MGEdge()) == MGVertex());
  BOOST_CHECK(G.GetSource(fail_edge) == MGVertex());
  BOOST_CHECK(G.GetTarget(MGEdge()) == MGVertex());
  BOOST_CHECK(G.GetTarget(fail_edge) == MGVertex());
  
  // Check that edge gets vertices correctly
  std::vector<std::pair<MGVertex, MGVertex>> expect_vertices, obtain_vertices;
  for (auto uv : incident_edges) expect_vertices.emplace_back(verts[uv.first], verts[uv.second]);
  for (auto e : edges) obtain_vertices.push_back(G.GetVertices(e));
  BOOST_CHECK_EQUAL_COLLECTIONS(expect_vertices.begin(), expect_vertices.end(),
                                obtain_vertices.begin(), obtain_vertices.end());
  
  // Check that edge gets vertices fails correctly
  BOOST_CHECK(G.GetVertices(MGEdge()) == std::make_pair(MGVertex(), MGVertex()));
  BOOST_CHECK(G.GetVertices(fail_edge) == std::make_pair(MGVertex(), MGVertex()));
}

BOOST_AUTO_TEST_CASE(clear) {
  G.Clear();
  // Check that vertices and edges are cleared
  BOOST_CHECK(G.NumEdges() == 0);
  BOOST_CHECK(G.NumVertices() == 0);
  auto vert_itr = G.GetVertices();
  auto edge_itr = G.GetEdges();
  BOOST_CHECK(vert_itr.first == vert_itr.second);
  BOOST_CHECK(edge_itr.first == edge_itr.second);
}

BOOST_AUTO_TEST_SUITE_END();

