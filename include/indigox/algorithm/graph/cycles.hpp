#include <algorithm>
#include <iterator>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <type_traits>
#include <vector>

#include "../../graph/base_graph.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP
#define INDIGOX_ALGORITHM_GRAPH_CYCLES_HPP

template <typename InputIter, typename T>
size_t combinations(InputIter first, InputIter last, size_t r,
                    std::vector<std::vector<T>>& out) {
  out.clear();
  std::vector<T> pool(first, last);
  if (r > pool.size()) return 0;
  std::vector<size_t> indices(r);
  std::iota(indices.begin(), indices.end(), 0);
  out.emplace_back(std::vector<T>());
  for (size_t i : indices) out.back().emplace_back(pool[i]);
  std::vector<size_t> tmp(r);
  std::iota(tmp.begin(), tmp.end(), 0);
  std::reverse(tmp.begin(), tmp.end());
  while (true) {
    size_t i = 0;
    bool broken = false;
    for (size_t x : tmp) {
      if (indices[x] != x + pool.size() - r) { i = x; broken = true; break; }
    }
    if (!broken) break;
    ++indices[i];
    std::vector<size_t> tmp2(r - i - 1);
    std::iota(tmp2.begin(), tmp2.end(), i + 1);
    for (size_t j : tmp2) indices[j] = indices[j - 1] + 1;
    out.emplace_back(std::vector<T>());
    for (size_t i : indices) out.back().emplace_back(pool[i]);
  }
  return out.size();
}

namespace indigox::algorithm {
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
  size_t CycleBasis(const graph::IXGraphBase<V, E, U, VP, EP>& G,
                   std::vector<std::vector<V*>>& paths) {
    paths.clear();
    std::vector<V*> verts;
    G.GetVertices(verts);
    while (!verts.empty()) {  // Loop over connected components
      V* root = verts.back(); verts.pop_back();
      std::vector<V*> stack; stack.push_back(root);
      std::map<V*, V*> pred; pred.emplace(root,root);
      std::map<V*, std::set<V*>> used; used.emplace(root, std::set<V*>());
      
      while (!stack.empty()) {  // walk the spanning tree finding cycles
        V* z = stack.back(); stack.pop_back();  // use last in so cycles easier
        std::set<V*>& zused = used.at(z);
        std::vector<V*> nbrs;
        G.GetNeighbours(z, nbrs);
        for (V* nbr : nbrs) {
          if (used.find(nbr) == used.end()) {   // new vert
            pred.emplace(nbr, z);
            stack.push_back(nbr);
            used.emplace(nbr, std::set<V*>{z});
          }
          // self loops
          else if (nbr == z) paths.emplace_back(std::vector<V*>{z});
          else if (zused.find(nbr) == zused.end()) {  // found a cycle
            std::set<V*>& pn = used.at(nbr);
            std::vector<V*> cycle{nbr, z};
            V* p = pred.at(z);
            while (pn.find(p) == pn.end()) {
              cycle.push_back(p);
              p = pred.at(p);
            }
            cycle.push_back(p);
            paths.emplace_back(cycle.begin(), cycle.end());
            used.at(nbr).insert(z);
          }
        }
      }
      auto p = [&pred] (V* v) { return pred.find(v) == pred.end(); };
      verts.erase(std::partition(verts.begin(), verts.end(), p), verts.end());
    }
    return paths.size();
  }
  
  template <class V, class E, class VP, class EP>
  size_t AllCycles(const graph::IXGraphBase<V, E, U, VP, EP>& G,
                 std::vector<std::vector<E*>>& cycles) {
    static_assert(!std::is_null_pointer<E>::value, "Edge type required");
    std::vector<std::vector<V*>> basis;
    std::vector<std::set<E*>> basis_edges;
    std::vector<E*> cyclic_edges;
    CycleBasis(G, basis);
    for (auto& path : basis) {
      path.push_back(path.front());
      std::set<E*> edges;
      for (size_t i = 1; i < path.size(); ++i)
        edges.insert(G.GetEdge(path[i-1], path[i]));
      std::vector<E*> tmp;
      std::set_union(edges.begin(), edges.end(),
                     cyclic_edges.begin(), cyclic_edges.end(),
                     std::back_inserter(tmp));
      cyclic_edges.swap(tmp);
      basis_edges.emplace_back(edges.begin(), edges.end());
    }
    std::vector<size_t> indices(basis_edges.size());
    std::iota(indices.begin(), indices.end(), 0);
    for (size_t r = 1; r < basis_edges.size() + 1; ++r) {
      std::vector<std::vector<size_t>> combos;
      combinations(indices.begin(), indices.end(), r, combos);
      for (auto& combo : combos) {
        std::vector<E*> xor_set;
        while (!combo.empty()) {
          size_t i = combo.back();
          combo.pop_back();
          std::vector<E*> tmp;
          std::set_symmetric_difference(xor_set.begin(), xor_set.end(),
                                        basis_edges[i].begin(),
                                        basis_edges[i].end(),
                                        std::back_inserter(tmp));
          xor_set.swap(tmp);
        }
        cycles.emplace_back(xor_set.begin(), xor_set.end());
      }
    }
    return cycles.size();
  }
  
  template <class V, class E, class VP, class EP>
  size_t MinimalCycles(const graph::IXGraphBase<V, E, U, VP, EP>& G,
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
