#include "../../utils/fwd_declares.hpp"

#include <map>
#include <vector>

#ifndef INDIGOX_ALGORITHM_GRAPH_PATHS_HPP
#define INDIGOX_ALGORITHM_GRAPH_PATHS_HPP

namespace indigox::algorithm {

  /*! \brief Find the shortest path between two vertices.
   *  \details Determines the shortest path between two vertices, using a
   *  bi-directional breadth-first search method. If there is no path between
   *  the vertices, the returned path is empty. The path is returned has the
   *  vertices to traverse in order from the source to the target. Multiple
   *  short paths may exist. This function returns only one of them.
   *  \tparam GraphType the type of the graph.
   *  \param G the graph to find the path in.
   *  \param source the source vertex.
   *  \param target the target vertex.
   *  \return the shortest path from \p source to \p target.
   *  \throws std::runtime_error if source and target are the same vertex or
   *  either of them is not part of the graph. */
  template <class V, class E, class S, class D, class VP, class EP>
  std::vector<E> ShortestPath(graph::BaseGraph<V, E, S, D, VP, EP> &G, V source,
                              V target);

  /*! \brief Find all the simple paths between two vertices.
   *  \details Determines all the paths between two vertices, using a modified
   *  depth-first search method. The order of the found paths is arbitary,
   *  though each path is ordered from source to target.
   *  \param G the graph to find paths in.
   *  \param source the source vertex.
   *  \param target the target vertex.
   *  \param[out] paths vector to store all the found paths in.
   *  \param limit the maximum length path to return.
   *  \throws std::runtime_error if you try something stupid. */
  template <class V, class E, class S, class D, class VP, class EP>
  void AllSimplePaths(graph::BaseGraph<V, E, S, D, VP, EP> &G, V source,
                      V target, std::vector<std::vector<E>> &paths,
                      int64_t limit = -1);

  template <class V> struct TraversalResults {
    using OrderType = typename std::vector<V>;
    using PredType = typename std::map<V, V>;
    using LengthType = typename std::map<V, int32_t>;
    OrderType discover_order;
    PredType predecessors;
    LengthType path_lengths;
    V furthest;

    TraversalResults(OrderType &order, PredType &pred, LengthType &length,
                     V &far)
        : discover_order(order), predecessors(pred), path_lengths(length),
          furthest(far) {}
  };

  /*! \brief Perform a breadth first search of G.
   *  \details Returns ordered vector of vertex pairs. Order is the order in
   *  which the vertices were discovered. First member of pair is the discovered
   *  vertex, second member is the predecessor of the discovered vertex. If no
   *  source is provided, the source is chosen arbitarily until all vertices
   *  have been discovered. If source is provided, only vertices within the
   *  component containing source will be searched. Limit is the limit of search
   *  depth.
   */
  template <class V, class E, class S, class D, class VP, class EP>
  TraversalResults<V> DepthFirstSearch(graph::BaseGraph<V, E, S, D, VP, EP> &G,
                                       V source = V(), int64_t limit = -1);

  template <class V, class E, class S, class D, class VP, class EP>
  TraversalResults<V>
  BreadthFirstSearch(graph::BaseGraph<V, E, S, D, VP, EP> &G, V source = V(),
                     int64_t limit = -1);
} // namespace indigox::algorithm

#endif /* INDIGOX_ALGORITHM_GRAPH_PATHS_HPP */
