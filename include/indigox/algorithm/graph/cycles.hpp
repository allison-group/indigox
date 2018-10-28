#include <algorithm>
#include <iterator>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <type_traits>
#include <vector>

#include "../../graph/base_graph.hpp"
#include "../../utils/numerics.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP
#define INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP

namespace indigox::algorithm {
  
  template<class V, class E, class D, class VP, class EP>
  int64_t CycleBasis(graph::BaseGraph<V,E,D,VP,EP>& G,
                     typename graph::BaseGraph<V,E,D,VP,EP>::CycleVertContain& basis) {
    static_assert(!D::is_directed, "Cycles require an undirected graph");
    using VertSet = eastl::vector_set<V>;
    using VertContain = typename graph::BaseGraph<V,E,D,VP,EP>::VertContain;
    
    basis.clear();
    VertContain verts = G.GetVertices();
    while (!verts.empty()) {
      V root = verts.back(); verts.pop_back();
      VertContain stack; stack.push_back(root);
      eastl::vector_map<V, V> pred; pred.emplace(root, root);
      eastl::vector_map<V, VertSet> used;
      used.emplace(root, VertSet());
      
      while (!stack.empty()) {  // walk the spanning tree finding cycles
        V z = stack.back(); stack.pop_back();  // use last in so cycles easier
        for (const V& nbr : G.GetNeighbours(z)) {
          if (used.find(nbr) == used.end()) {   // new vert
            if (pred.find(nbr) == pred.end()) pred.emplace(nbr, z);
            else pred.at(nbr) = z;
            stack.push_back(nbr);
            used.emplace(nbr, VertSet{z});
          }
          // self loops
          else if (nbr == z) basis.emplace_back(VertContain{z});
          else if (used.at(z).find(nbr) == used.at(z).end()) {  // found a cycle
            VertSet& pn = used.at(nbr);
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
      auto part = [&pred] (const V& v) { return pred.find(v) == pred.end(); };
      verts.erase(std::partition(verts.begin(), verts.end(), part), verts.end());
    }
    return basis.size();
  }
  
  template<class V, class E, class D, class VP, class EP>
  int64_t CycleBasis(graph::BaseGraph<V,E,D,VP,EP>& G,
                     typename graph::BaseGraph<V,E,D,VP,EP>::CycleEdgeContain& e_basis) {
    using CycleVertContain = typename graph::BaseGraph<V,E,D,VP,EP>::CycleVertContain;
    using VertContain = typename graph::BaseGraph<V,E,D,VP,EP>::VertContain;
    using EdgeContain = typename graph::BaseGraph<V,E,D,VP,EP>::EdgeContain;
    
    CycleVertContain v_basis;
    CycleBasis(G, v_basis);
    e_basis.clear();
    if (v_basis.size() == 0) return 0;
    for (VertContain& path : v_basis) {
      path.emplace(path.end(), path.back());
      EdgeContain path_edge;
      for (size_t i = 1; i < path.size(); ++i)
        path_edge.emplace(path_edge.end(), G.GetEdge(path[i-1], path[i]));
      e_basis.emplace(e_basis.end(), path_edge.begin(), path_edge.end());
    }
    return e_basis.size();
  }
  
  template<class V, class E, class D, class VP, class EP>
  int64_t AllCycles(graph::BaseGraph<V,E,D,VP,EP>& G,
                    typename graph::BaseGraph<V,E,D,VP,EP>::CycleEdgeContain& ecycles) {
    
    using CycleEdgeContain = typename graph::BaseGraph<V,E,D,VP,EP>::CycleEdgeContain;
    using EdgeContain = typename graph::BaseGraph<V,E,D,VP,EP>::EdgeContain;
    
    CycleEdgeContain basis_edge;
    CycleBasis(G, basis_edge);
    std::vector<size_t> indices(basis_edge.size());
    std::iota(indices.begin(), indices.end(), 0);
    for (size_t r = 1; r <= basis_edge.size(); ++r) {
      std::vector<std::vector<size_t>> combs;
      Combinations(indices.begin(), indices.end(), r, combs);
      for (std::vector<size_t>& combo : combs) {
        EdgeContain xor_set, tmp;
        while (!combo.empty()) {
          tmp.clear();
          size_t i = combo.back(); combo.pop_back();
          std::set_symmetric_difference(xor_set.begin(), xor_set.end(),
                                        basis_edge[i].begin(),
                                        basis_edge[i].end(),
                                        std::back_inserter(tmp));
          xor_set.swap(tmp);
        }
        ecycles.emplace(ecycles.end(), xor_set.begin(), xor_set.end());
      }
    }
    return ecycles.size();
  }
  
  using U = graph::Undirected;
  
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
  
  template <class V, class E, class VP, class EP>
  size_t MinimalCycles(const graph::BaseGraph<V, E, U, VP, EP>& G,
                      std::vector<std::vector<E*>>& min_cycles,
                      bool strict = false) {
    // Gets all the cycles then removes those it doesn't need.
    std::vector<std::vector<E*>> all_cycles;
    AllCycles(G, all_cycles);
    std::set<E*> cyclic_edges;
    for (auto& cycle : all_cycles) {
      std::set<E*> tmp;
      std::set_union(cyclic_edges.begin(), cyclic_edges.end(),
                     cycle.begin(), cycle.end(),
                     std::inserter(tmp, tmp.begin()));
      cyclic_edges.swap(tmp);
    }
    
    // Sort so smallest cycles are first
    auto sorter = [](const std::vector<E*>& a, const std::vector<E*>& b) {
      return a.size() < b.size();
    };
    std::stable_sort(all_cycles.begin(), all_cycles.end(), sorter);
    
    size_t current_size = all_cycles.front().size();
    std::map<E*, size_t> edge_counts;
    for (auto& cycle : all_cycles) {
      /*  Only check for all edges being accounted for when the cycles is larger
       *  than the previous one, unless the strict flag is set. */
      if (strict || cycle.size() > current_size) {
        std::set<E*> added_edges;
        for (auto& c : min_cycles) {
          std::set<E*> tmp;
          std::set_union(added_edges.begin(), added_edges.end(),
                         c.begin(), c.end(),
                         std::inserter(tmp, tmp.begin()));
          added_edges.swap(tmp);
        }
        if (added_edges == cyclic_edges) break;
        current_size = cycle.size();
      }
      
      // Don't add cycles which will mean more than 2 occurances of an edge
      bool add_cycle = true;
      for (E* e : cycle) {
        if (edge_counts.find(e) == edge_counts.end()) continue;
        if (edge_counts[e] < 2) continue;
        add_cycle = false;
        break;
      }
      
      // Don't add cycles that don't add new edges to the total set
      if (add_cycle) {
        bool new_edge = false;
        for (E* e : cycle) {
          if (edge_counts.find(e) != edge_counts.end()) continue;
          new_edge = true;
          break;
        }
        if (!new_edge) add_cycle = false;
      }
      
      if (add_cycle) {
        for (E* e : cycle) {
          if (edge_counts.find(e) == edge_counts.end())
            edge_counts.emplace(e,1);
          else ++edge_counts[e];
        }
        min_cycles.emplace_back(cycle.begin(), cycle.end());
      }
    }
    return min_cycles.size();
  }
  
}

#endif /* INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP */
