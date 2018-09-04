#include <map>
#include <numeric>
#include <queue>
#include <vector>

#include <EASTL/vector_map.h>
#include <EASTL/vector_set.h>

#include "../../utils/fwd_declares.hpp"
#include "../../utils/numerics.hpp"

#ifndef INDIGOX_ALGORITHM_GRAPH_PATHS_HPP
#define INDIGOX_ALGORITHM_GRAPH_PATHS_HPP

namespace indigox::algorithm {
  
  template <class VertType>
  struct Path : public std::vector<VertType> {
    using std::vector<VertType>::vector;
    typename std::vector<VertType>::size_type Length() const {
      return std::vector<VertType>::size() ? std::vector<VertType>::size() - 1: 0;
    }
  };
  
  template <class EdgeType>
  struct EdgePath : public std::vector<EdgeType> {
    using std::vector<EdgeType>::vector;
    typename std::vector<EdgeType>::size_type Length() const {
      return std::vector<EdgeType>::size();
    }
  };
  
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
  template <class GraphType>
  Path<typename GraphType::VertexType>
  ShortestPath(const std::shared_ptr<GraphType>& G,
               const typename GraphType::VertexType& source,
               const typename GraphType::VertexType& target) {
    using Vertex = typename GraphType::VertexType;
    using PathType = Path<Vertex>;
    using MapType = eastl::vector_map<Vertex, Vertex>;
    
    if (!G) throw std::runtime_error("Graph is empty");
    if (!G->HasVertex(source)) throw std::runtime_error("Source not in graph");
    if (!G->HasVertex(target)) throw std::runtime_error("Target not in graph");
    if (source == target) throw std::runtime_error("Source and target are the same");
    
    MapType predecessors, successors;
    predecessors.emplace(source, Vertex());
    successors.emplace(target, Vertex());
    
    std::vector<Vertex> forward, backward;
    forward.push_back(source);
    backward.push_back(target);
    Vertex midpoint = source;
    bool path_found = false;
    
    // Find the path
    while (!forward.empty() && !backward.empty()) {
      if (path_found) break;
      std::vector<Vertex> current;
      if (forward.size() <= backward.size()) {
        current.swap(forward);
        for (Vertex v : current) {
          for (auto nbrs = G->GetNeighbours(v); nbrs.first != nbrs.second; ++nbrs.first) {
            Vertex n = *nbrs.first;
            if (predecessors.find(n) == predecessors.end()) {
              forward.push_back(n);
              predecessors.emplace(n, v);
            }
            if (successors.find(n) != successors.end()) {
              midpoint = n;
              path_found = true;
              break;
            }
          }
          if (path_found) break;
        }
      } else {
        current.swap(backward);
        for (Vertex v : current) {
          for (auto nbrs = G->GetNeighbours(v); nbrs.first != nbrs.second; ++nbrs.first) {
            Vertex n = *nbrs.first;
            if (successors.find(n) == successors.end()) {
              backward.push_back(n);
              successors.emplace(n, v);
            }
            if (predecessors.find(n) != predecessors.end()) {
              midpoint = n;
              path_found = true;
              break;
            }
          }
          if (path_found) break;
        }
      }
    }
    if (!path_found) return PathType();
    
    // Build the path from source to midpoint
    PathType path; path.reserve(predecessors.size() + successors.size());
    while (midpoint) {
      path.push_back(midpoint);
      midpoint = predecessors.at(midpoint);
    }
    std::reverse(path.begin(), path.end());
    midpoint = successors.at(path.back());
    while (midpoint) {
      path.push_back(midpoint);
      midpoint = successors.at(midpoint);
    }
    return path;
  }
  
  
  /*! \brief Find the shortest path between two vertices.
   *  \details Determines the shortest path between two vertices, using a
   *  bi-directional breadth-first search method. If there is no path between
   *  the vertices, the returned path is empty. The path is returned has the
   *  edges to traverse in order from the source to the target. Multiple short
   *  paths may exist. This function returns only one of them.
   *  \tparam GraphType the type of the graph.
   *  \param G the graph to find the path in.
   *  \param source the source vertex.
   *  \param target the target vertex.
   *  \return the shortest path from \p source to \p target.
   *  \throws std::runtime_error if source and target are the same vertex or
   *  either of them is not part of the graph. */
  template <class GraphType>
  EdgePath<typename GraphType::EdgeType>
  ShortestEdgePath(const std::shared_ptr<GraphType>& G,
                   const typename GraphType::VertexType& source,
                   const typename GraphType::VertexType& target) {
    using Vertex = typename GraphType::VertexType;
    using Edge = typename GraphType::EdgeType;
    using PathType = Path<Vertex>;
    using EdgePathType = EdgePath<Edge>;
    
    PathType p = ShortestPath(G, source, target);
    if (p.size() < 2) return EdgePathType();
    EdgePathType path; path.reserve(p.size());
    for (size_ i = 0; i < p.size() - 1; ++i)
      path.push_back(G->GetEdge(p[i], p[i+1]));
    return path;
  }
  
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
  template <class GraphType>
  void AllSimplePaths(const std::shared_ptr<GraphType>& G,
                      const typename GraphType::VertexType& source,
                      const typename GraphType::VertexType& target,
                      std::vector<Path<typename GraphType::VertexType>>& paths,
                      size_ limit = std::numeric_limits<size_>::max()) {
    using Vertex = typename GraphType::VertexType;
    using NbrsIters = typename GraphType::NbrsIter;
    using NbrsPair = std::pair<NbrsIters, NbrsIters>;
    
    if (!G) throw std::runtime_error("Graph is empty");
    if (!G->HasVertex(source)) throw std::runtime_error("Source not in graph");
    if (!G->HasVertex(target)) throw std::runtime_error("Target not in graph");
    if (source == target) throw std::runtime_error("Source and target are the same");
    
    paths.clear();
    
    std::vector<Vertex> vis;
    vis.reserve(G->NumVertices());
    vis.push_back(source);
    
    std::vector<NbrsPair> stack;
    stack.reserve(G->NumVertices());
    stack.emplace_back(G->GetNeighbours(source));
    
    while (!stack.empty()) {
      NbrsPair& cs = stack.back();
      if (cs.first == cs.second) { // if child is None
        stack.pop_back();
        vis.pop_back();
        continue;
      }
      Vertex c = *cs.first;
      ++cs.first;
      if (vis.size() < limit) { // elif len(visited < cutoff
        if (c == target) {
          paths.emplace_back(vis.begin(), vis.end());
          paths.back().emplace_back(target);
        } else if (std::find(vis.begin(), vis.end(), c) == vis.end()) {
          vis.push_back(c);
          stack.emplace_back(G->GetNeighbours(c));
        }
      } else {
        if (c == target || std::find(cs.first, cs.second, target) != cs.second) {
          paths.emplace_back(vis.begin(), vis.end());
          paths.back().emplace_back(target);
        }
        stack.pop_back();
        vis.pop_back();
      }
    }
  }
  
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
  template <class GraphType>
  void AllSimpleEdgePaths(const std::shared_ptr<GraphType>& G,
                          const typename GraphType::VertexType& source,
                          const typename GraphType::VertexType& target,
                          std::vector<EdgePath<typename GraphType::EdgeType>>& paths,
                          size_ limit = std::numeric_limits<size_>::max()) {
    using Vertex = typename GraphType::VertexType;
    using Edge = typename GraphType::EdgeType;
    using PathType = Path<Vertex>;
    using EdgePathType = EdgePath<Edge>;
    
    std::vector<PathType> p;
    AllSimplePaths(G, source, target, p, limit);
    for (PathType& path : p) {
      paths.emplace_back(EdgePathType());
      for (size_ i = 0; i < path.size() - 1; ++i)
        paths.back().push_back(G->GetEdge(path[i], path[i+1]));
    }
  }
  
 
}

#endif /* INDIGOX_ALGORITHM_GRAPH_PATHS_HPP */
