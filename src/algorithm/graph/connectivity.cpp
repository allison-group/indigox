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
#include <memory>
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

  // =======================================================================
  // == Connected subgraphs implementation =================================
  // =======================================================================

  template <class GraphType> struct ConnectedSubgraphs<GraphType>::Impl {

    using vert_contain = typename GraphType::VertContain;
    using BitSet = boost::dynamic_bitset<>;
    using StackItem = stdx::triple<BitSet>;

    GraphType graph;
    size_t min_subgraph_size;
    size_t max_subgraph_size;
    vert_contain vertices;
    eastl::vector_map<size_t, BitSet> neighbours;
    std::vector<StackItem> stack;

    Impl(GraphType &G, size_t min, size_t max)
        : graph(G), min_subgraph_size(min), max_subgraph_size(max),
          vertices(G.GetVertices()) {
    }

    void BuildNeighboursBitsets() {
      for (size_t i = 0; i < vertices.size(); ++i) {
        auto &v_nbrs = graph.GetNeighbours(vertices[i]);
        BitSet nbrs(vertices.size());
        nbrs.reset();
        for (auto &v : v_nbrs) {
          auto pos = std::find(vertices.begin(), vertices.end(), v);
          nbrs.set(std::distance(vertices.begin(), pos));
        }
        neighbours.emplace(i, nbrs);
      }
    }

    void Initialise() {
      BuildNeighboursBitsets();
      BitSet bag(vertices.size());
      bag.reset();
      for (size_t i = 0; i < vertices.size(); ++i)
        bag.set(i);
      BitSet initial(vertices.size());
      initial.reset();
      BitSet nbrs(bag);
      nbrs.reset();
      size_t pos = initial.find_first();
      while (pos < neighbours.size()) {
        auto test = neighbours.find(pos);
        if (test == neighbours.end())
          break;
        nbrs |= test->second;
        pos = initial.find_next(pos);
      }
      stack.emplace_back(bag, initial, nbrs);
    }

    bool NextSubgraph(GraphType &subgraph) {
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

        if (possible.none() && cur_subg.any() &&
            cur_subg.count() <= max_subgraph_size &&
            cur_subg.count() >= min_subgraph_size) {
          std::vector<typename GraphType::VertexType> verts;
          verts.reserve(cur_subg.count());
          size_t pos = cur_subg.find_first();
          while (pos < vertices.size()) {
            verts.emplace_back(vertices[pos]);
            pos = cur_subg.find_next(pos);
          }
          GraphType subg = graph.Subgraph(verts);
          subgraph = subg;
          return true;
        } else if (possible.any()) {
          BitSet v = possible;
          v.reset();
          v.set(possible.find_first());
          if (cur_subg.count() <= max_subgraph_size) {
            BitSet bag_minus_v = cur_bag;
            bag_minus_v.reset(v.find_first());
            stack.emplace_back(bag_minus_v, cur_subg, cur_nbrs);
            stack.emplace_back(bag_minus_v, cur_subg | v,
                               cur_nbrs | neighbours.at(v.find_first()));
          }
        }
      }
      return false;
    }
  };

  template <class GraphType>
  ConnectedSubgraphs<GraphType>::ConnectedSubgraphs(GraphType &G, size_t min,
                                                    size_t max)
      : implementation(std::make_unique<Impl>(G, min, max)) {
    implementation->Initialise();
  }

  template <class GraphType>
  ConnectedSubgraphs<GraphType>::~ConnectedSubgraphs() = default;

  template <class GraphType>
  bool ConnectedSubgraphs<GraphType>::operator()(GraphType &subgraph) {
    return implementation->NextSubgraph(subgraph);
  }

  template class ConnectedSubgraphs<graph::MolecularGraph>;
  template class ConnectedSubgraphs<graph::CondensedMolecularGraph>;

} // namespace indigox::algorithm
