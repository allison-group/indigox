#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iterator>
#include <random>
#include <set>
#include <vector>

#include "class_test_wrappers.hpp"

namespace indigox::test {
  struct AssignmentGraphFixture {
    Molecule mol, rand_mol;
    indigox::test::IXAssignmentGraph g, g_rand;
    AssignmentGraphFixture() : mol(CreateMolecule()), rand_mol(CreateMolecule()),
    g(mol->GetGraph()), g_rand(rand_mol->GetGraph()) {
      // create the simple molecule
      Atom c = mol->NewAtom(); Atom h1 = mol->NewAtom(); Atom h2 = mol->NewAtom();
      Atom h3 = mol->NewAtom(); Atom h4 = mol->NewAtom();
      mol->NewBond(c, h1); mol->NewBond(c, h2); mol->NewBond(c, h3);
      mol->NewBond(c, h4);
      g = indigox::test::IXAssignmentGraph(mol->GetGraph());
      
      // create a random molecule
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<size_> at_dis(15,35);
      std::uniform_int_distribution<size_> bn_dis(35,55);
      
      size_ num = at_dis(gen);
      std::vector<Atom> atoms;
      mol->ReserveAtoms(num);
      atoms.reserve(num);
      for (size_ i = 0; i < num; ++i) atoms.emplace_back(rand_mol->NewAtom());
      
      std::vector<std::pair<size_, size_>> possible_bonds, samples;
      possible_bonds.reserve(num * (num + 1) / 2);
      for (size_ i = 0; i < num; ++i) {
        for (size_ j = i + 1; j < num; ++j ) possible_bonds.emplace_back(i,j);
      }

      num = bn_dis(gen);
      std::sample(possible_bonds.begin(), possible_bonds.end(),
                  std::back_inserter(samples), num, gen);
      for (auto& i : samples) rand_mol->NewBond(atoms[i.first], atoms[i.second]);
      g_rand = indigox::test::IXAssignmentGraph(rand_mol->GetGraph());
      
    }
  };
}

using namespace indigox;
using namespace indigox::graph;

BOOST_FIXTURE_TEST_SUITE(ixassignment_graph, indigox::test::AssignmentGraphFixture);

BOOST_AUTO_TEST_CASE(constructor) {
  BOOST_CHECK(g.NumVertices() == (mol->NumAtoms() + mol->NumBonds()));
  BOOST_CHECK(g_rand.NumVertices() == (rand_mol->NumAtoms() + rand_mol->NumBonds()));
  auto g_verts = g.GetVertices();
  auto ran_verts = g_rand.GetVertices();
  BOOST_CHECK(g_verts.second == (g_verts.first + g.NumVertices()));
  BOOST_CHECK(ran_verts.second == (ran_verts.first + g_rand.NumVertices()));
}

BOOST_AUTO_TEST_CASE(has_vertex) {
  // Check that all rand_mol vert/edges correctly have
  std::vector<bool> expected(rand_mol->NumBonds() + rand_mol->NumAtoms(), true);
  std::vector<bool> obtained;
  obtained.reserve(expected.size());
  MolecularGraph molG = rand_mol->GetGraph();
  for (auto it = molG->GetVertices(); it.first != it.second; ++it.first)
    obtained.emplace_back(g_rand.HasVertex(*it.first));
  for (auto it = molG->GetEdges(); it.first != it.second; ++it.first)
    obtained.emplace_back(g_rand.HasVertex(*it.first));
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                obtained.begin(), obtained.end());
  
  // Check a newly added vertex is correctly haved
  MGVertex v = molG->GetVertex(rand_mol->NewAtom());
  AGVertex u = g_rand.AddVertex(v);
  BOOST_CHECK(g_rand.HasVertex(u));
  
  // Check that correct failing
  // fail on null passes
  BOOST_CHECK(!g_rand.HasVertex(AGVertex()));
  BOOST_CHECK(!g_rand.HasVertex(MGVertex()));
  BOOST_CHECK(!g_rand.HasVertex(MGEdge()));
  // fail on not part of passes
  BOOST_CHECK(!g_rand.HasVertex(*g.GetVertices().first));
  BOOST_CHECK(!g_rand.HasVertex(*mol->GetGraph()->GetVertices().first));
  BOOST_CHECK(!g_rand.HasVertex(*mol->GetGraph()->GetEdges().first));
}

BOOST_AUTO_TEST_CASE(get_vertex) {
  // check correct gets
  auto g_verts = g_rand.GetVertices();
  std::set<AGVertex> expected(g_verts.first, g_verts.second);
  std::set<AGVertex> obtained;
  auto m_verts = rand_mol->GetGraph()->GetVertices();
  auto m_edges = rand_mol->GetGraph()->GetEdges();
  for (; m_verts.first != m_verts.second; ++m_verts.first)
    obtained.emplace(g_rand.GetVertex(*m_verts.first));
  for (; m_edges.first != m_edges.second; ++m_edges.first)
    obtained.emplace(g_rand.GetVertex(*m_edges.first));

  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                obtained.begin(), obtained.end());
  
  // check correct failing
  // fail on null passes
  BOOST_CHECK(!g_rand.GetVertex(MGVertex()));
  BOOST_CHECK(!g_rand.GetVertex(MGEdge()));
  // fail on not part of passes
  BOOST_CHECK(!g_rand.GetVertex(*mol->GetGraph()->GetVertices().first));
  BOOST_CHECK(!g_rand.GetVertex(*mol->GetGraph()->GetEdges().first));
}

BOOST_AUTO_TEST_CASE(degree) {
  // check degrees are correct
  MolecularGraph m = rand_mol->GetGraph();
  std::vector<size_> expected, obtained;
  expected.reserve(m->NumVertices() + m->NumEdges());
  obtained.reserve(m->NumVertices() + m->NumEdges());
  for (auto it = g_rand.GetVertices(); it.first != it.second; ++it.first) {
    AGVertex v = *it.first;
    if (v->IsVertexMapped()) {
      expected.emplace_back(m->Degree(v->GetSourceVertex()));
      obtained.emplace_back(g_rand.Degree(v));
    } else {
      expected.emplace_back(2);
      obtained.emplace_back(g_rand.Degree(v));
    }
  }
  
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                obtained.begin(), obtained.end());
  
  // check correct failing
  size_ num = std::numeric_limits<size_>::max();
  BOOST_CHECK(g_rand.Degree(AGVertex()) == num);  // fail on null passes
  BOOST_CHECK(g_rand.Degree(*g.GetVertices().first) == num);  // fail on not part
  
}

BOOST_AUTO_TEST_CASE(connected) {
  BOOST_CHECK(g.IsConnected());
}

BOOST_AUTO_TEST_CASE(vertex_assigned) {
  std::vector<uint_> expected, obtained;
  expected.reserve(g_rand.NumVertices() * 2);
  obtained.reserve(g_rand.NumVertices() * 2);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint_> dis(0,10);
  for (auto it = g_rand.GetVertices(); it.first != it.second; ++it.first) {
    AGVertex v = *it.first;
    uint_ pre = dis(gen);
    uint_ post = dis(gen);
    v->SetPreAssignedCount(pre);
    v->SetAssignedCount(post);
    expected.emplace_back(pre);
    expected.emplace_back(pre + post);
    obtained.emplace_back(v->GetPreAssignedCount());
    obtained.emplace_back(v->GetTotalAssigned());
  }
  
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                obtained.begin(), obtained.end());
  
  // check correctly inits to 0
  expected.clear();
  expected.assign(g.NumVertices() * 2, 0);
  obtained.clear();
  for (auto it = g.GetVertices(); it.first != it.second; ++it.first) {
    obtained.emplace_back((*it.first)->GetPreAssignedCount());
    obtained.emplace_back((*it.first)->GetTotalAssigned());
  }
  
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                obtained.begin(), obtained.end());
}

BOOST_AUTO_TEST_CASE(vertex_sources) {
  auto m_verts = rand_mol->GetGraph()->GetVertices();
  auto m_edges = rand_mol->GetGraph()->GetEdges();
  std::set<MGVertex> expected_vert(m_verts.first, m_verts.second);
  std::set<MGEdge> expected_edge(m_edges.first, m_edges.second);
  std::set<MGEdge> obtained_edge;
  std::set<MGVertex> obtained_vert;
  for (auto it = g_rand.GetVertices(); it.first != it.second; ++it.first) {
    AGVertex v = *it.first;
    if (v->IsEdgeMapped()) {
      obtained_edge.emplace(v->GetSourceEdge());
      BOOST_CHECK(!v->GetSourceVertex());
    }
    else if (v->IsVertexMapped()) {
      obtained_vert.emplace(v->GetSourceVertex());
      BOOST_CHECK(!v->GetSourceEdge());
    }
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_vert.begin(), expected_vert.end(),
                                obtained_vert.begin(), obtained_vert.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_edge.begin(), expected_edge.end(),
                                obtained_edge.begin(), obtained_edge.end());
}

BOOST_AUTO_TEST_CASE(neighbours) {
  // Check all neighbours are correct
  MolecularGraph m = rand_mol->GetGraph();
  for (auto it = g_rand.GetVertices(); it.first != it.second; ++it.first) {
    std::set<AGVertex> expected, obtained;
    AGVertex v = *it.first;
    if (v->IsVertexMapped()) {
      MGVertex mv = v->GetSourceVertex();
      auto i = m->GetNeighbours(mv);
      std::vector<MGVertex> v_nbrs(i.first, i.second);
      for (MGVertex u : v_nbrs) {
        MGEdge e = m->GetEdge(u, mv);
        expected.emplace(g_rand.GetVertex(e));
      }
    } else {
      MGEdge e = v->GetSourceEdge();
      expected.emplace(g_rand.GetVertex(m->GetSource(e)));
      expected.emplace(g_rand.GetVertex(m->GetTarget(e)));
    }
    auto j = g_rand.GetNeighbours(v);
    obtained.insert(j.first, j.second);
    
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  obtained.begin(), obtained.end());
  }
  
  // Check correctly no neighbours when bad vertex
  auto it = g_rand.GetNeighbours(AGVertex()); // null pass
  BOOST_CHECK(std::distance(it.first, it.second) == 0);
  it = g_rand.GetNeighbours(*g.GetVertices().first);  // not part of pass
  BOOST_CHECK(std::distance(it.first, it.second) == 0);
}

BOOST_AUTO_TEST_SUITE_END();

