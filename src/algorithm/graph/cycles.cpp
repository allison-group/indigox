#include <indigox/algorithm/graph/cycles.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

#include <algorithm>
#include <iterator>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <type_traits>
#include <vector>

#include <EASTL/vector_set.h>

namespace indigox::algorithm {

  using namespace indigox::graph;

  // ===========================================================================
  // == Cycle vertex basis implementation ======================================
  // ===========================================================================

  template <class V, class E, class S, class D, class VP, class EP>
  int64_t
  VCycleBasis(BaseGraph<V, E, S, D, VP, EP> &G,
              typename BaseGraph<V, E, S, D, VP, EP>::CycleVertContain &basis) {
    static_assert(!D::is_directed, "Cycles require an undirected graph");
    using VertSet = eastl::vector_set<V>;
    using VertContain = typename BaseGraph<V, E, S, D, VP, EP>::VertContain;

    basis.clear();
    VertContain verts = G.GetVertices();
    while (!verts.empty()) {
      V root = verts.back();
      verts.pop_back();
      VertContain stack;
      stack.push_back(root);
      eastl::vector_map<V, V> pred;
      pred.emplace(root, root);
      eastl::vector_map<V, VertSet> used;
      used.emplace(root, VertSet());

      while (!stack.empty()) { // walk the spanning tree finding cycles
        V z = stack.back();
        stack.pop_back(); // use last in so cycles easier
        for (const V &nbr : G.GetNeighbours(z)) {
          if (used.find(nbr) == used.end()) { // new vert
            if (pred.find(nbr) == pred.end())
              pred.emplace(nbr, z);
            else
              pred.at(nbr) = z;
            stack.push_back(nbr);
            used.emplace(nbr, VertSet{z});
          }
          // self loops
          else if (nbr == z)
            basis.emplace_back(VertContain{z});
          else if (used.at(z).find(nbr) == used.at(z).end()) { // found a cycle
            VertSet &pn = used.at(nbr);
            VertContain cycle{nbr, z};
            V p = pred.at(z);
            while (pn.find(p) == pn.end()) {
              cycle.push_back(p);
              p = pred.at(p);
            }
            cycle.push_back(p);
            basis.emplace_back(cycle.begin(), cycle.end());
            used.at(nbr).insert(z);
          }
        }
      }
      auto part = [&pred](const V &v) { return pred.find(v) == pred.end(); };
      verts.erase(std::partition(verts.begin(), verts.end(), part),
                  verts.end());
    }
    return basis.size();
  }

  template <>
  int64_t CycleBasis(CondensedMolecularGraph::graph_type &G,
                     CondensedMolecularGraph::CycleVertContain &basis) {
    return VCycleBasis(G, basis);
  }

  template <>
  int64_t CycleBasis(MolecularGraph::graph_type &G,
                     MolecularGraph::CycleVertContain &basis) {
    return VCycleBasis(G, basis);
  }

  // ===========================================================================
  // == Cycle edge basis implementation ========================================
  // ===========================================================================

  template <class V, class E, class S, class D, class VP, class EP>
  int64_t ECycleBasis(
      BaseGraph<V, E, S, D, VP, EP> &G,
      typename BaseGraph<V, E, S, D, VP, EP>::CycleEdgeContain &e_basis) {
    using CycleVertContain =
        typename BaseGraph<V, E, S, D, VP, EP>::CycleVertContain;
    using VertContain = typename BaseGraph<V, E, S, D, VP, EP>::VertContain;
    using EdgeContain = typename BaseGraph<V, E, S, D, VP, EP>::EdgeContain;

    CycleVertContain v_basis;
    CycleBasis(G, v_basis);
    e_basis.clear();
    if (v_basis.size() == 0) return 0;
    for (VertContain &path : v_basis) {
      path.emplace(path.end(), path.front());
      EdgeContain path_edge;
      for (size_t i = 1; i < path.size(); ++i)
        path_edge.emplace(path_edge.end(), G.GetEdge(path[i - 1], path[i]));
      e_basis.emplace(e_basis.end(), path_edge.begin(), path_edge.end());
    }
    return e_basis.size();
  }

  template <>
  int64_t CycleBasis(CondensedMolecularGraph::graph_type &G,
                     CondensedMolecularGraph::CycleEdgeContain &basis) {
    return ECycleBasis(G, basis);
  }

  template <>
  int64_t CycleBasis(MolecularGraph::graph_type &G,
                     MolecularGraph::CycleEdgeContain &basis) {
    return ECycleBasis(G, basis);
  }

  // ===========================================================================
  // == All cycles implementation ==============================================
  // ===========================================================================

  template <class V, class E, class S, class D, class VP, class EP>
  bool __ring_connectivity_check(BaseGraph<V, E, S, D, VP, EP> &G,
                                 std::vector<E>& edges) {
    eastl::vector_set<V> all_vertices;
    all_vertices.reserve(edges.size());
    for (E e: edges) {
      all_vertices.insert(G.GetSourceVertex(e));
      all_vertices.insert(G.GetTargetVertex(e));
    }
    
    std::vector<V> tmp(all_vertices.begin(), all_vertices.end());
    S subG = G.Subgraph(tmp, edges);
    return subG.IsConnected();
  }
  
  template <class V, class E, class S, class D, class VP, class EP,
            class Container>
  int64_t AllCycles(BaseGraph<V, E, S, D, VP, EP> &G, Container &ecycles) {
    using CycleEdgeContain =
        typename BaseGraph<V, E, S, D, VP, EP>::CycleEdgeContain;
    using EdgeContain = typename BaseGraph<V, E, S, D, VP, EP>::EdgeContain;

    ecycles.clear();
    CycleEdgeContain basis_edge;
    CycleBasis(G, basis_edge);
    
    //  Group indices into connected groups.
    //  Only do combinations of connected groups
    std::vector<size_t> indices(basis_edge.size());
    std::iota(indices.begin(), indices.end(), 0);
    eastl::vector_set<size_t> idxs(indices.begin(), indices.end());
    
    std::vector<std::vector<size_t>> index_grouping;
    while (idxs.size()) {
      index_grouping.push_back(std::vector<size_t>());
      index_grouping.back().push_back(idxs.back()); idxs.pop_back();
      bool added_new = true;
      while (added_new) {
        added_new = false;
        std::vector<E> current_e;
        for (size_t i : index_grouping.back()) current_e.insert(current_e.end(), basis_edge[i].begin(), basis_edge[i].end());
        for (size_t i : idxs) {
          std::vector<E> test_e(current_e.begin(), current_e.end());
          test_e.insert(test_e.end(), basis_edge[i].begin(), basis_edge[i].end());
          if (__ring_connectivity_check(G, test_e)) {
            added_new = true;
            idxs.erase(i);
            index_grouping.back().push_back(i);
            break;
          }
        }
      }
    }
    
    for (std::vector<size_t>& group : index_grouping) {
      for (size_t r = 1; r <= group.size(); ++r) {
        std::vector<std::vector<size_t>> combs;
        Combinations(group.begin(), group.end(), r, combs);
        for (std::vector<size_t> &combo : combs) {
          EdgeContain xor_set, tmp;
          while (!combo.empty()) {
            tmp.clear();
            size_t i = combo.back();
            combo.pop_back();
            std::set_symmetric_difference(
                                          xor_set.begin(), xor_set.end(), basis_edge[i].begin(),
                                          basis_edge[i].end(), std::back_inserter(tmp));
            xor_set.swap(tmp);
          }
          
          // Only add the cycle if its connected
          if (__ring_connectivity_check(G, xor_set))
            ecycles.emplace(ecycles.end(), xor_set.begin(), xor_set.end());
        }
      }
    }
    return ecycles.size();
  }

  template int64_t AllCycles(CondensedMolecularGraph::graph_type &,
                             CondensedMolecularGraph::CycleEdgeContain &);
  template int64_t AllCycles(MolecularGraph::graph_type &,
                             MolecularGraph::CycleEdgeContain &);

  /*! \brief Calculate a list of cycles which form a basis for cycles of G.
   *  \details A basis for cycles of a graph is a minimal collection of cycles
   *  such that any cycle in the graph can be written as a sum of cycles in the
   *  basis, where a summation of cycles is defined as an XOR of the edges. This
   *  is only implemented for Undirected graphs.
   *
   *  This implementation follows that of the networkx python library.
   *  \tparam V, E, VP, EP template parameters of the graph.
   *  \param G the graph to determine cycles of.
   *  \param[out] paths vector of paths defining each member of cycle basis.
   */

  template <class V, class E, class S, class U, class VP, class EP>
  size_t MinimalCycles(const graph::BaseGraph<V, E, S, U, VP, EP> &G,
                       std::vector<std::vector<E *>> &min_cycles,
                       bool strict = false) {
    // Gets all the cycles then removes those it doesn't need.
    std::vector<std::vector<E *>> all_cycles;
    AllCycles(G, all_cycles);
    std::set<E *> cyclic_edges;
    for (auto &cycle : all_cycles) {
      std::set<E *> tmp;
      std::set_union(cyclic_edges.begin(), cyclic_edges.end(), cycle.begin(),
                     cycle.end(), std::inserter(tmp, tmp.begin()));
      cyclic_edges.swap(tmp);
    }

    // Sort so smallest cycles are first
    auto sorter = [](const std::vector<E *> &a, const std::vector<E *> &b) {
      return a.size() < b.size();
    };
    std::stable_sort(all_cycles.begin(), all_cycles.end(), sorter);

    size_t current_size = all_cycles.front().size();
    std::map<E *, size_t> edge_counts;
    for (auto &cycle : all_cycles) {
      /*  Only check for all edges being accounted for when the cycles is larger
       *  than the previous one, unless the strict flag is set. */
      if (strict || cycle.size() > current_size) {
        std::set<E *> added_edges;
        for (auto &c : min_cycles) {
          std::set<E *> tmp;
          std::set_union(added_edges.begin(), added_edges.end(), c.begin(),
                         c.end(), std::inserter(tmp, tmp.begin()));
          added_edges.swap(tmp);
        }
        if (added_edges == cyclic_edges) break;
        current_size = cycle.size();
      }

      // Don't add cycles which will mean more than 2 occurances of an edge
      bool add_cycle = true;
      for (E *e : cycle) {
        if (edge_counts.find(e) == edge_counts.end()) continue;
        if (edge_counts[e] < 2) continue;
        add_cycle = false;
        break;
      }

      // Don't add cycles that don't add new edges to the total set
      if (add_cycle) {
        bool new_edge = false;
        for (E *e : cycle) {
          if (edge_counts.find(e) != edge_counts.end()) continue;
          new_edge = true;
          break;
        }
        if (!new_edge) add_cycle = false;
      }

      if (add_cycle) {
        for (E *e : cycle) {
          if (edge_counts.find(e) == edge_counts.end())
            edge_counts.emplace(e, 1);
          else
            ++edge_counts[e];
        }
        min_cycles.emplace_back(cycle.begin(), cycle.end());
      }
    }
    return min_cycles.size();
  }

} // namespace indigox::algorithm
