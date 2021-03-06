#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

#include <EASTL/vector_map.h>
#include <deque>
#include <vector>

namespace indigox::algorithm {

  using namespace indigox::graph;

  // ===========================================================================
  // == ShortestPath implementation ============================================
  // ===========================================================================
  template <class V, class E, class S, class D, class VP, class EP>
  std::vector<E> ShortestPath(BaseGraph<V, E, S, D, VP, EP> &G, V source,
                              V target) {
    if (!G.HasVertex(source) || !G.HasVertex(target))
      throw std::runtime_error("Vertices not in graph");
    std::vector<E> path;
    if (source == target) return path;

    using NbrType = eastl::vector_map<V, V>;
    NbrType pre, suc;
    pre.emplace(source, V());
    suc.emplace(target, V());
    std::vector<V> forward, backward;
    forward.emplace_back(source);
    backward.emplace_back(target);
    V midpoint = source;
    bool found = false;

    // Find the vertices to traverse in the path
    while (!forward.empty() && !backward.empty()) {
      if (found) break;
      std::vector<V> current_dir;
      if (forward.size() <= backward.size()) {
        current_dir.swap(forward);
        for (V v : current_dir) {
          if (found) break;
          for (V nbr : G.GetNeighbours(v)) {
            if (pre.find(nbr) == pre.end()) {
              forward.emplace_back(nbr);
              pre.emplace(nbr, v);
            }
            if (suc.find(nbr) != suc.end()) {
              midpoint = nbr;
              found = true;
              break;
            }
          }
        }
      } else {
        current_dir.swap(backward);
        for (V v : current_dir) {
          if (found) break;
          for (V nbr : G.GetNeighbours(v)) {
            if (suc.find(nbr) == suc.end()) {
              backward.emplace_back(nbr);
              suc.emplace(nbr, v);
            }
            if (pre.find(nbr) != pre.end()) {
              midpoint = nbr;
              found = true;
              break;
            }
          }
        }
      }
    }
    if (!found) return path;

    // Build the path
    path.reserve(pre.size() + suc.size());
    V mid = midpoint;
    while (pre.at(mid)) {
      V previous = pre.at(mid);
      path.emplace_back(G.GetEdge(previous, mid));
      mid = previous;
    }
    std::reverse(path.begin(), path.end());
    while (suc.at(midpoint)) {
      V next = suc.at(midpoint);
      path.emplace_back(G.GetEdge(midpoint, next));
      midpoint = next;
    }
    return path;
  }

  template std::vector<CMGEdge>
  ShortestPath(CondensedMolecularGraph::graph_type &, CMGVertex, CMGVertex);

  template std::vector<MGEdge> ShortestPath(MolecularGraph::graph_type &,
                                            MGVertex, MGVertex);

  // ===========================================================================
  // == AllSimplePaths implementation ==========================================
  // ===========================================================================

  template <class V, class E, class S, class D, class VP, class EP>
  void AllSimplePaths(BaseGraph<V, E, S, D, VP, EP> &G, V source, V target,
                      std::vector<std::vector<E>> &paths, int64_t limit) {
    using NbrsType = typename BaseGraph<V, E, S, D, VP, EP>::VertContain;
    using PathType = std::vector<E>;
    paths.clear();
    if (!G.HasVertex(source) || !G.HasVertex(target))
      throw std::runtime_error("Vertices not in graph");
    if (source == target) return;

    std::vector<V> visited;
    visited.reserve(G.NumVertices());
    visited.push_back(source);

    std::vector<NbrsType> stack;
    stack.reserve(G.NumVertices());
    stack.emplace_back(G.GetNeighbours(source));

    while (!stack.empty()) {
      NbrsType &children = stack.back();
      if (children.empty()) {
        stack.pop_back();
        visited.pop_back();
        continue;
      }

      V child = children.back();
      children.pop_back();
      if (limit < 0 || static_cast<int64_t>(visited.size()) < limit) {
        if (child == target) {
          paths.emplace_back(PathType());
          for (size_t i = 1; i < visited.size(); ++i)
            paths.back().emplace_back(G.GetEdge(visited[i - 1], visited[i]));
          paths.back().emplace_back(G.GetEdge(visited.back(), target));
        } else if (std::find(visited.begin(), visited.end(), child) ==
                   visited.end()) {
          visited.emplace_back(child);
          stack.emplace_back(G.GetNeighbours(child));
        }
      } else {
        if (child == target || std::find(children.begin(), children.end(),
                                         target) != children.end()) {
          paths.emplace_back(PathType());
          for (size_t i = 1; i < visited.size(); ++i)
            paths.back().emplace_back(G.GetEdge(visited[i - 1], visited[i]));
          paths.back().emplace_back(G.GetEdge(visited.back(), target));
        }
        stack.pop_back();
        visited.pop_back();
      }
    }
  }

  template void AllSimplePaths(CondensedMolecularGraph::graph_type &, CMGVertex,
                               CMGVertex, std::vector<std::vector<CMGEdge>> &,
                               int64_t);
  template void AllSimplePaths(MolecularGraph::graph_type &, MGVertex, MGVertex,
                               std::vector<std::vector<MGEdge>> &, int64_t);

  // ===========================================================================
  // == DepthFirstSearch implementation ========================================
  // ===========================================================================

  template <class V, class E, class S, class D, class VP, class EP>
  TraversalResults<V> DepthFirstSearch(graph::BaseGraph<V, E, S, D, VP, EP> &G,
                                       V source, int64_t limit) {
    using GraphType = typename graph::BaseGraph<V, E, S, D, VP, EP>;
    using NbrsIter =
        typename GraphType::NbrsContain::mapped_type::const_iterator;
    using Results = TraversalResults<V>;

    std::vector<V> vertices;
    if (!source)
      vertices.insert(vertices.end(), G.GetVertices().begin(),
                      G.GetVertices().end());
    else if (G.HasVertex(source))
      vertices.emplace_back(source);
    else
      throw std::runtime_error("Source vertex not part of graph");

    eastl::vector_set<V> visited;
    typename Results::OrderType discover_order;
    typename Results::PredType predecessors;
    typename Results::LengthType path_lengths;
    V furthest;

    discover_order.reserve(G.NumVertices());

    if (limit < 1) limit = G.NumVertices();
    for (V start : vertices) {
      if (visited.find(start) != visited.end()) continue;
      visited.insert(start);
      discover_order.emplace_back(start);
      predecessors[start] = V();
      path_lengths[start] = 0;
      furthest = start;
      using StackItem = stdx::triple<V, int64_t, std::pair<NbrsIter, NbrsIter>>;
      std::vector<StackItem> stack;
      stack.emplace_back(start, limit,
                         std::make_pair(G.GetNeighbours(start).begin(),
                                        G.GetNeighbours(start).end()));
      while (!stack.empty()) {
        StackItem &item = stack.back();
        if (item.third.first == item.third.second) {
          stack.pop_back();
          continue;
        }
        V child = *item.third.first;
        ++item.third.first;
        if (visited.find(child) != visited.end()) continue;
        discover_order.emplace_back(child);
        predecessors[child] = item.first;
        path_lengths[child] = 1 + path_lengths[item.first];
        if (path_lengths[child] > path_lengths[furthest]) furthest = child;
        visited.insert(child);
        if (item.second > 1)
          stack.emplace_back(child, item.second - 1,
                             std::make_pair(G.GetNeighbours(child).begin(),
                                            G.GetNeighbours(child).end()));
      }
    }

    return Results(discover_order, predecessors, path_lengths, furthest);
  }

  template TraversalResults<CMGVertex>
  DepthFirstSearch(CondensedMolecularGraph::graph_type &, CMGVertex, int64_t);
  template TraversalResults<MGVertex>
  DepthFirstSearch(MolecularGraph::graph_type &, MGVertex, int64_t);

  // ===========================================================================
  // == BreadthFirstSearch implementation ======================================
  // ===========================================================================

  template <class V, class E, class S, class D, class VP, class EP>
  TraversalResults<V>
  BreadthFirstSearch(graph::BaseGraph<V, E, S, D, VP, EP> &G, V source,
                     int64_t limit) {
    using GraphType = typename graph::BaseGraph<V, E, S, D, VP, EP>;
    using NbrsIter =
        typename GraphType::NbrsContain::mapped_type::const_iterator;
    using Results = TraversalResults<V>;

    std::vector<V> vertices;
    if (!source)
      vertices.insert(vertices.end(), G.GetVertices().begin(),
                      G.GetVertices().end());
    else if (G.HasVertex(source))
      vertices.emplace_back(source);
    else
      throw std::runtime_error("Source vertex not part of graph");

    eastl::vector_set<V> visited;
    typename Results::OrderType discover_order;
    typename Results::PredType predecessors;
    typename Results::LengthType path_lengths;
    V furthest;

    discover_order.reserve(G.NumVertices());

    if (limit < 1) limit = G.NumVertices();

    for (V start : vertices) {
      if (visited.find(start) != visited.end()) continue;
      visited.insert(start);
      discover_order.emplace_back(start);
      predecessors[start] = V();
      path_lengths[start] = 0;
      furthest = start;

      using QueueItem = stdx::triple<V, int64_t, std::pair<NbrsIter, NbrsIter>>;
      std::deque<QueueItem> queue;
      queue.emplace_back(start, limit,
                         std::make_pair(G.GetNeighbours(start).begin(),
                                        G.GetNeighbours(start).end()));
      while (!queue.empty()) {
        QueueItem &item = queue.front();
        if (item.third.first == item.third.second) {
          queue.pop_front();
          continue;
        }
        V child = *item.third.first;
        ++item.third.first;
        if (visited.find(child) != visited.end()) continue;
        discover_order.emplace_back(child);
        predecessors[child] = item.first;
        path_lengths[child] = 1 + path_lengths[item.first];
        if (path_lengths[child] > path_lengths[furthest]) furthest = child;
        visited.insert(child);
        if (item.second > 1)
          queue.emplace_back(child, item.second - 1,
                             std::make_pair(G.GetNeighbours(child).begin(),
                                            G.GetNeighbours(child).end()));
      }
    }

    return Results(discover_order, predecessors, path_lengths, furthest);
  }

  template TraversalResults<CMGVertex>
  BreadthFirstSearch(CondensedMolecularGraph::graph_type &, CMGVertex, int64_t);

  template TraversalResults<MGVertex>
  BreadthFirstSearch(MolecularGraph::graph_type &, MGVertex, int64_t);
} // namespace indigox::algorithm
