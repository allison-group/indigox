#include <indigox/algorithm/access.hpp>
#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/triple.hpp>

#include <boost/dynamic_bitset.hpp>
#include <boost/graph/connected_components.hpp>

#include <EASTL/bitset.h>
#include <EASTL/vector_map.h>
#include <vector>

namespace indigox::algorithm {

  using namespace indigox::graph;

  // ===========================================================================
  // == Connected components implementation ====================================
  // ===========================================================================

  template <class V, class E, class S, class D, class VP, class EP,
            class Container>
  int64_t ConnectedComponents(BaseGraph<V, E, S, D, VP, EP> &G,
                              Container &bits) {
    static_assert(!D::is_directed, "Connectivity requires an undirected graph");

    using GraphType = BaseGraph<V, E, S, D, VP, EP>;
    using BaseType = typename GraphType::graph_type;
    using VertIdxMap = eastl::vector_map<typename GraphType::VertType, int>;
    using namespace boost;

    BaseType &g = access::GetGraph(G);
    VertIdxMap data;
    associative_property_map<VertIdxMap> indexMap(data);
    typename GraphType::VertIter begin, end;
    std::tie(begin, end) = boost::vertices(g);
    for (int i = 0; begin != end; ++begin, ++i)
      put(indexMap, *begin, i);
    size_t num = connected_components(g, get(&VP::component, g),
                                      vertex_index_map(indexMap));
    bits.clear();
    bits.assign(num, typename GraphType::VertContain());
    for (auto &v : access::GetVertexMap(G).right)
      bits[g[v.first].component].emplace_back(v.second);
    return bits.size();
  }

  template int64_t
  ConnectedComponents(CondensedMolecularGraph::graph_type &,
                      CondensedMolecularGraph::ComponentContain &);
  template int64_t ConnectedComponents(MolecularGraph::graph_type &,
                                       MolecularGraph::ComponentContain &);

  // ===========================================================================
  // == Connected subgraphs implementation =====================================
  // ===========================================================================

  template <class BitSet>
  void __bitset_connected_subgraphs(BitSet bag, BitSet initial_subg,
                                    eastl::vector_map<size_t, BitSet> &all_nbrs,
                                    std::vector<BitSet> &subGs, size_t minsize,
                                    size_t maxsize) {
    subGs.clear();

    // Determine the current neighbours of a bag
    BitSet nbrs(bag);
    nbrs.reset();
    size_t pos = initial_subg.find_first();
    while (pos < all_nbrs.size()) {
      if (all_nbrs.find(pos) == all_nbrs.end())
        break;
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
      if (cur_subg.none())
        possible = cur_bag;
      else
        possible = cur_bag & cur_nbrs;

      if (possible.none() && cur_subg.any() && cur_subg.count() <= maxsize &&
          cur_subg.count() >= minsize)
        subGs.emplace_back(cur_subg);
      else if (possible.any()) {
        BitSet v = possible;
        v.reset();
        v.set(possible.find_first());
        if (cur_subg.count() <= maxsize) {
          BitSet bag_minus_v = cur_bag;
          bag_minus_v.reset(v.find_first());
          stack.emplace_back(bag_minus_v, cur_subg, cur_nbrs);
          stack.emplace_back(bag_minus_v, cur_subg | v,
                             cur_nbrs | all_nbrs.at(v.find_first()));
        }
      }
    }
  }

  template <class BitSet, class V, class E, class S, class D, class VP,
            class EP>
  void __build_neighbours_bitsets(
      graph::BaseGraph<V, E, S, D, VP, EP> &G,
      const typename graph::BaseGraph<V, E, S, D, VP, EP>::VertContain
          &vertices,
      eastl::vector_map<size_t, BitSet> &nbrs) {
    for (size_t i = 0; i < vertices.size(); ++i) {
      auto &v_nbrs = G.GetNeighbours(vertices[i]);
      BitSet n(vertices.size());
      n.reset();
      for (auto &v : v_nbrs) {
        auto pos = std::find(vertices.begin(), vertices.end(), v);
        n.set(std::distance(vertices.begin(), pos));
      }
      nbrs.emplace(i, n);
    }
  }

  template <class BitSet, class V, class E, class S, class D, class VP,
            class EP>
  void __populate_subgraphs(
      BaseGraph<V, E, S, D, VP, EP> &G,
      const typename BaseGraph<V, E, S, D, VP, EP>::VertContain &vertices,
      std::vector<BitSet> &subg, std::vector<S> &subGs) {
    for (BitSet &sub : subg) {
      std::vector<V> verts;
      verts.reserve(vertices.size());
      size_t pos = sub.find_first();
      while (pos < vertices.size()) {
        verts.emplace_back(vertices[pos]);
        pos = sub.find_next(pos);
      }
      subGs.emplace_back(G.Subgraph(verts));
    }
  }

  template <class BitSet, class V, class E, class S, class D, class VP,
            class EP>
  void _runner(BaseGraph<V, E, S, D, VP, EP> &G, std::vector<S> &sub_graphs,
               size_t minsize, size_t maxsize) {
    using GraphType = BaseGraph<V, E, S, D, VP, EP>;
    using VertContain = typename GraphType::VertContain;
    sub_graphs.clear();

    const VertContain &vertices = G.GetVertices();
    // Build up the neighbours
    eastl::vector_map<size_t, BitSet> nbrs;
    __build_neighbours_bitsets(G, vertices, nbrs);

    // Calculate the subgs
    std::vector<BitSet> subg;
    BitSet bag(vertices.size());
    bag.reset();
    for (size_t i = 0; i < vertices.size(); ++i)
      bag.set(i);
    BitSet initial(vertices.size());
    initial.reset();
    __bitset_connected_subgraphs(bag, initial, nbrs, subg, minsize, maxsize);

    // Populate the subGs
    __populate_subgraphs(G, vertices, subg, sub_graphs);
  }

  template <class V, class E, class S, class D, class VP, class EP>
  int64_t ConnectedSubgraphs(BaseGraph<V, E, S, D, VP, EP> &G,
                             std::vector<S> &subs, size_t min, size_t max) {
    using Bitset32 = eastl::bitset<32, uint32_t>;
    using Bitset64 = eastl::bitset<64, uint64_t>;
    using Bitset128 = eastl::bitset<128, uint64_t>;
    using BitsetMax = boost::dynamic_bitset<>;
    if (G.NumVertices() <= 32)
      _runner<Bitset32>(G, subs, min, max);
    else if (G.NumVertices() <= 64)
      _runner<Bitset64>(G, subs, min, max);
    else if (G.NumVertices() <= 128)
      _runner<Bitset128>(G, subs, min, max);
    else
      _runner<BitsetMax>(G, subs, min, max);
    return subs.size();
  }

  template int64_t ConnectedSubgraphs(CondensedMolecularGraph::graph_type &,
                                      std::vector<CondensedMolecularGraph> &,
                                      size_t, size_t);
  template int64_t ConnectedSubgraphs(MolecularGraph::graph_type &,
                                      std::vector<MolecularGraph> &, size_t,
                                      size_t);

} // namespace indigox::algorithm
