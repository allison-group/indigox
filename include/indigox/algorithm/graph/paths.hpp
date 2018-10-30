#include <vector>

#include "../../utils/fwd_declares.hpp"

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
  std::vector<E> ShortestPath(graph::BaseGraph<V,E,S,D,VP,EP>& G,
                              V source, V target);
  
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
  void AllSimplePaths(graph::BaseGraph<V,E,S,D,VP,EP>& G, V source, V target,
                      std::vector<std::vector<E>>& paths, int64_t limit = -1);
}

#endif /* INDIGOX_ALGORITHM_GRAPH_PATHS_HPP */
