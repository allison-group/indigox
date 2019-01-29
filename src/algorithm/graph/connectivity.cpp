#include <indigox/algorithm/access.hpp>
#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/molecule.hpp>
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

  // =======================================================================
  // == Optimal charge groups implementation ===============================
  // =======================================================================
  struct ChargeGroupOptimiser {
    using BitSet = boost::dynamic_bitset<>;
    using Score = double;
    using Division = std::vector<BitSet>;
    using ScoredDivisions = std::map<Division, Score>;
    using OptimalDivision = std::map<BitSet, Division>;

    graph::MolecularGraph full_graph, leafless;
    int32_t size_limit;

    std::vector<graph::MGVertex> non_leaf_vertices, leaf_vertices;
    std::vector<int32_t> weights;

    std::map<BitSet::size_type, Score> position_scores;
    std::map<BitSet::size_type, BitSet> neighbours;

    ScoredDivisions seen_division_scores;
    OptimalDivision optimal_seen_divisions;

    ChargeGroupOptimiser(const Molecule &mol, int32_t limit);

    double PenaltyScore(BitSet group);

    void Initalise();

    std::vector<std::vector<graph::MGVertex>> Optimise();

    Division &Divisions(BitSet bag);
    void ConnectedSubgraphs(BitSet bag, BitSet current, BitSet nbrs,
                            std::vector<BitSet> &subgraphs);
  };

  std::vector<std::vector<Atom>> OptimalChargeGroups(const Molecule &mol,
                                                     int32_t limit) {
    ChargeGroupOptimiser optimiser(mol, limit);
    optimiser.Initalise();
    std::vector<std::vector<graph::MGVertex>> v_cgs = optimiser.Optimise();
    std::vector<std::vector<Atom>> charge_groups;
    for (std::vector<graph::MGVertex> grp : v_cgs) {
      charge_groups.emplace_back(std::vector<Atom>());
      for (graph::MGVertex v : grp) {
        charge_groups.back().emplace_back(v.GetAtom());
      }
    }
    return charge_groups;
  }

  ChargeGroupOptimiser::ChargeGroupOptimiser(const Molecule &mol, int32_t limit)
      : full_graph(mol.GetGraph()), size_limit(limit) {
  }

  double ChargeGroupOptimiser::PenaltyScore(BitSet group) {
    double score = 0.;
    BitSet::size_type pos = group.find_first();
    while (pos < group.size()) {
      score += position_scores[pos];
      pos = group.find_next(pos);
    }
    return abs(score);
  }

  void ChargeGroupOptimiser::Initalise() {
    // Vertices identify and sort
    for (graph::MGVertex v : full_graph.GetVertices()) {
      if (full_graph.Degree(v) == 1) {
        leaf_vertices.push_back(v);
      } else {
        non_leaf_vertices.push_back(v);
      }
    }
    leafless = full_graph.Subgraph(non_leaf_vertices);

    BitSet::size_type v_pos = 0;
    for (graph::MGVertex v : non_leaf_vertices) {
      // Neighbours generation

      BitSet nbrs(non_leaf_vertices.size());
      for (graph::MGVertex n : leafless.GetNeighbours(v)) {
        BitSet::size_type n_pos = std::distance(
            non_leaf_vertices.begin(),
            std::find(non_leaf_vertices.begin(), non_leaf_vertices.end(), n));
        nbrs.set(n_pos);
      }
      neighbours.emplace(v_pos, nbrs);

      // Scores per vertex pre calculation
      int32_t weight = 1;
      Atom atm = v.GetAtom();
      double score = atm.GetFormalCharge() - atm.GetPartialCharge();
      for (graph::MGVertex u : full_graph.GetNeighbours(v)) {
        if (std::find(non_leaf_vertices.begin(), non_leaf_vertices.end(), u) ==
            non_leaf_vertices.end()) {
          atm = u.GetAtom();
          score += atm.GetFormalCharge() - atm.GetPartialCharge();
          ++weight;
        }
      }
      position_scores.emplace(v_pos, score);
      weights.emplace_back(weight);
      ++v_pos;

      // Terminal case seen score
      seen_division_scores.emplace(Division(), 0.0);
      optimal_seen_divisions.emplace(BitSet(non_leaf_vertices.size()),
                                     Division());
    }
  }

  std::vector<std::vector<graph::MGVertex>> ChargeGroupOptimiser::Optimise() {
    BitSet bag(non_leaf_vertices.size());
    bag.set();
    Division optimal_division = Divisions(bag);

    std::vector<std::vector<graph::MGVertex>> charge_groups;

    for (BitSet div : optimal_division) {
      charge_groups.emplace_back(std::vector<graph::MGVertex>());
      charge_groups.back().reserve(div.count());
      BitSet::size_type pos = div.find_first();
      while (pos < div.size()) {
        charge_groups.back().emplace_back(non_leaf_vertices[pos]);
        for (graph::MGVertex nbr :
             full_graph.GetNeighbours(non_leaf_vertices[pos])) {
          if (full_graph.Degree(nbr) == 1) {
            charge_groups.back().emplace_back(nbr);
          }
        }
        pos = div.find_next(pos);
      }
    }

    return charge_groups;
  }

  ChargeGroupOptimiser::Division &ChargeGroupOptimiser::Divisions(BitSet bag) {
    auto seen_pos = optimal_seen_divisions.find(bag);
    if (seen_pos != optimal_seen_divisions.end()) {
      return seen_pos->second;
    }

    BitSet::size_type v = bag.find_first();
    std::vector<BitSet> subgraphs;
    BitSet bag_minus_v(bag);
    bag_minus_v.reset(v);
    BitSet current(bag.size()), nbrs(neighbours.at(v));
    current.set(v);
    ConnectedSubgraphs(bag_minus_v, current, nbrs, subgraphs);

    for (BitSet g : subgraphs) {
      BitSet bag_minus_g = bag - g;
      Division new_div = Divisions(bag_minus_g);
      Score new_score = seen_division_scores.at(new_div);
      Score cur_score = PenaltyScore(g);

      auto seen_position = optimal_seen_divisions.find(bag);
      Division optimal_div(new_div);
      optimal_div.emplace_back(g);
      Score optimal_score = cur_score + new_score;
      if (seen_position == optimal_seen_divisions.end()) {
        seen_division_scores.emplace(optimal_div, optimal_score);
        optimal_seen_divisions.emplace(bag, optimal_div);
      } else if ((seen_division_scores.at(seen_position->second) >
                  (optimal_score))) {
        seen_position->second = optimal_div;
        seen_division_scores[optimal_div] = optimal_score;
      }

      if (-1e-10 < optimal_score && 1e-10 > optimal_score)
        break;
    }
    return optimal_seen_divisions.at(bag);
  }

  void ChargeGroupOptimiser::ConnectedSubgraphs(
      BitSet bag, BitSet current, BitSet nbrs, std::vector<BitSet> &subgraphs) {
    subgraphs.clear();
    using StackItem = stdx::triple<BitSet>;
    std::vector<StackItem> stack;
    stack.emplace_back(bag, current, nbrs);
    while (!stack.empty()) {
      StackItem item = stack.back();
      stack.pop_back();

      BitSet possibles = item.first;
      if (item.second.any())
        possibles = item.first & item.third;
      if (possibles.none() && item.second.any()) {
        subgraphs.emplace_back(item.second);
      } else if (possibles.any()) {
        BitSet::size_type v = possibles.find_first();
        int32_t v_weight = 0;
        BitSet::size_type pos = item.second.find_first();
        while (pos < item.second.size()) {
          v_weight += weights[pos];
          pos = item.second.find_next(pos);
        }
        if (v_weight < size_limit) {
          BitSet bag_minus_v(item.first);
          bag_minus_v.reset(v);
          BitSet g_and_v(item.second);
          g_and_v.set(v);
          BitSet new_nbrs = item.third | neighbours.at(v);
          stack.emplace_back(bag_minus_v, item.second, item.third);
          stack.emplace_back(bag_minus_v, g_and_v, new_nbrs);
        } else {
          subgraphs.emplace_back(item.second);
        }
      }
    }
    std::sort(subgraphs.begin(), subgraphs.end(),
              [](BitSet &a, BitSet &b) { return a.count() < b.count(); });
  }

} // namespace indigox::algorithm
