#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <EASTL/bitset.h>
#include <EASTL/vector_map.h>

#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/utils/triple.hpp>

namespace indigox::algorithm {
  
  template <class BitSet>
  void __bitset_connected_subgraphs(BitSet bag, BitSet initial_subg,
                                    eastl::vector_map<size_t, BitSet>& all_nbrs,
                                    std::vector<BitSet>& subGs,
                                    size_t minsize, size_t maxsize) {
    subGs.clear();
    
    // Determine the current neighbours of a bag
    BitSet nbrs(bag); nbrs.reset();
    size_t pos = initial_subg.find_first();
    while (pos < all_nbrs.size()) {
      if (all_nbrs.find(pos) == all_nbrs.end()) break;
      nbrs |= all_nbrs.at(pos);
      pos = initial_subg.find_next(pos);
    }
    
    // Make a stack
    using StackItem = stdx::triple<BitSet, BitSet, BitSet>;
    std::vector<StackItem> stack;
    stack.emplace_back(bag, initial_subg, nbrs);
    
    // Run through the stack
    while (stack.size()) {
      StackItem item = stack.back();
      stack.pop_back();
      BitSet cur_bag = item.first;
      BitSet cur_subg = item.second;
      BitSet cur_nbrs = item.third;
      
      BitSet possible;
      if (cur_subg.none()) possible = cur_bag;
      else possible = cur_bag & cur_nbrs;
      
      if (possible.none() && cur_subg.any() && cur_subg.count() <= maxsize
          && cur_subg.count() >= minsize) subGs.emplace_back(cur_subg);
      else if (possible.any()) {
        BitSet v = possible; v.reset();
        v.set(possible.find_first());
        if (cur_subg.count() <= maxsize) {
          BitSet bag_minus_v = cur_bag; bag_minus_v.reset(v.find_first());
          stack.emplace_back(bag_minus_v, cur_subg, cur_nbrs);
          stack.emplace_back(bag_minus_v, cur_subg | v,
                             cur_nbrs | all_nbrs.at(v.find_first()));
        }
      }
    }
  }
  
  template <class BitSet, class GraphType, class VertexType>
  void __build_neighbours_bitsets(const GraphType& G,
                                  const std::vector<VertexType>& vertices,
                                  eastl::vector_map<size_t, BitSet>& nbrs) {
    for (size_t i = 0; i < vertices.size(); ++i) {
      auto v_nbrs = G->GetNeighbours(vertices[i]);
      BitSet n(vertices.size()); n.reset();
      for (; v_nbrs.first != v_nbrs.second; ++v_nbrs.first) {
        auto pos = std::find(vertices.begin(), vertices.end(), *v_nbrs.first);
        n.set(std::distance(vertices.begin(), pos));
      }
      nbrs.emplace(i, n);
    }
  }
  
  template <class BitSet, class GraphType, class VertexType>
  void __populate_subgraphs(const GraphType& G,
                            const std::vector<VertexType>& vertices,
                            std::vector<BitSet>& subg,
                            std::vector<GraphType>& subGs) {
    for (BitSet& sub : subg) {
      std::vector<VertexType> verts;
      verts.reserve(vertices.size());
      size_t pos = sub.find_first();
      while (pos < vertices.size()) {
        verts.emplace_back(vertices[pos]);
        pos = sub.find_next(pos);
      }
      subGs.emplace_back(G->InduceSubgraph(verts.begin(), verts.end()));
    }
  }
  
  template <class BitSet, class GraphType, class VertexType>
  void __connected_subgraphs_runner(const GraphType& G,
                                    std::vector<GraphType>& subGs,
                                    size_t minsize, size_t maxsize) {
    subGs.clear();
    std::vector<VertexType> vertices(G->GetVertices().first,
                                     G->GetVertices().second);
    // Build up the neighbours
    eastl::vector_map<size_t, BitSet> nbrs;
    __build_neighbours_bitsets(G, vertices, nbrs);
    // Calculate the subgs
    std::vector<BitSet> subg;
    BitSet bag(vertices.size()); bag.reset();
    for (size_t i = 0; i < vertices.size(); ++i) bag.set(i);
    BitSet initial(vertices.size()); initial.reset();
    __bitset_connected_subgraphs(bag, initial, nbrs, subg, minsize, maxsize);
    // Populate the subGs
    __populate_subgraphs(G, vertices, subg, subGs);
  }
  
  size_t ConnectedSubgraphs(const graph::CondensedMolecularGraph& G,
                           std::vector<graph::CondensedMolecularGraph>& subGs,
                           size_t minsize, size_t maxsize) {
    if (G->NumVertices() <= 32)
      __connected_subgraphs_runner<eastl::bitset<32, uint32_t>,
      graph::CondensedMolecularGraph,
      graph::CMGVertex>(G, subGs, minsize, maxsize);
    else if (G->NumVertices() <= 64)
      __connected_subgraphs_runner<eastl::bitset<64, uint64_t>,
      graph::CondensedMolecularGraph,
      graph::CMGVertex>(G, subGs, minsize, maxsize);
    else
      __connected_subgraphs_runner<boost::dynamic_bitset<>,
      graph::CondensedMolecularGraph,
      graph::CMGVertex>(G, subGs, minsize, maxsize);
    
    return subGs.size();
  }
  
  size_t ConnectedSubgraphs(const graph::MolecularGraph& G,
                           std::vector<graph::MolecularGraph>& subGs,
                           size_t minsize, size_t maxsize) {
    if (G->NumVertices() <= 32)
      __connected_subgraphs_runner<eastl::bitset<32, uint32_t>,
                                   graph::MolecularGraph,
                                   graph::MGVertex>(G, subGs, minsize, maxsize);
    else if (G->NumVertices() <= 64)
      __connected_subgraphs_runner<eastl::bitset<64, uint64_t>,
                                   graph::MolecularGraph,
                                   graph::MGVertex>(G, subGs, minsize, maxsize);
    else if (G->NumVertices() <= 128)
      __connected_subgraphs_runner<eastl::bitset<128, uint64_t>,
                                   graph::MolecularGraph,
                                   graph::MGVertex>(G, subGs, minsize, maxsize);
    else
      __connected_subgraphs_runner<boost::dynamic_bitset<>,
                                   graph::MolecularGraph,
                                   graph::MGVertex>(G, subGs, minsize, maxsize);
    
    return subGs.size();
  }
}
