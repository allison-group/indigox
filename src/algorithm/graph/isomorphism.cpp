#include <boost/graph/vf2_sub_graph_iso.hpp>

#include <EASTL/vector_map.h>

#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

namespace indigox::algorithm {
  using namespace indigox::graph;
  
  template<class GraphType, class UserCallback>
  struct InternalCallback {
    using Graph = std::shared_ptr<GraphType>;
    using Vertex = typename GraphType::VertexType;
    using BaseGraph = typename GraphType::graph_type;
    using BoostGraph = typename BaseGraph::graph_type;
    using BoostVertex = typename BoostGraph::vertex_descriptor;
    using BoostEdge = typename BoostGraph::edge_descriptor;
    using BoostVertexIdxMap = eastl::vector_map<BoostVertex, size_>;
  public:
    Graph G_small, G_large;
    BaseGraph &base_small, &base_large;
    BoostVertexIdxMap &vmap;
    UserCallback& user_callback;
    
    InternalCallback(Graph small, Graph large, BoostVertexIdxMap& idxmap,
                     UserCallback& callback)
    : G_small(small), G_large(large), base_small(access::member_graph(small)),
    base_large(access::member_graph(large)), vmap(idxmap),
    user_callback(callback) { }
    
    template <class CMap1to2, class CMap2to1>
    bool operator()(CMap1to2 one, CMap2to1) {
      eastl::vector_map<Vertex, Vertex> map;
      typename BoostGraph::vertex_iterator b, e;
      std::tie(b, e) = boost::vertices(access::member_graph(base_small));
      for (; b != e; ++b) {
        auto vsmall = base_small.GetVertex(*b);
        auto vlarge = base_large.GetVertex(boost::get(one, *b));
        map.emplace(vsmall->shared_from_this(), vlarge->shared_from_this());
      }
      return user_callback(map);
    }
    
    bool operator()(BoostVertex, BoostVertex) { return true; }
    bool operator()(BoostEdge, BoostEdge) { return true; }
  };
  
  template <class UserCallback>
  struct CMGCallback
  : public InternalCallback<IXCondensedMolecularGraph, UserCallback> {
    using BaseType = InternalCallback<IXCondensedMolecularGraph, UserCallback>;
    using Graph = typename BaseType::Graph;
    using BaseGraph = typename BaseType::BaseGraph;
    using BoostGraph = typename BaseType::BoostGraph;
    using BoostVertex = typename BaseType::BoostVertex;
    using BoostEdge = typename BaseType::BoostEdge;
    using BoostVertexIdxMap = typename BaseType::BoostVertexIdxMap;
    
    using VertexIso = graph::VertexIsoMask;
    using EdgeIso = graph::EdgeIsoMask;
    using VertexMasks = eastl::vector_map<BoostVertex, VertexIso>;
    using EdgeMasks = eastl::vector_map<BoostEdge, EdgeIso>;
    
    VertexMasks vmasks;
    EdgeMasks emasks;
    
    CMGCallback(Graph small, Graph large, BoostVertexIdxMap& idxmap,
                UserCallback& callback)
    : BaseType(small, large, idxmap, callback) {
      graph::IXCondensedMolecularGraph::VertIter b, e;
      for (std::tie(b, e) = small->GetVertices(); b != e; ++b) {
        BoostVertex v = this->base_small.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetIsomorphismMask());
      }
      for (std::tie(b, e) = large->GetVertices(); b != e; ++b) {
        BoostVertex v = this->base_large.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetIsomorphismMask());
      }
      graph::IXCondensedMolecularGraph::EdgeIter B, E;
      for (std::tie(B, E) = small->GetEdges(); B != E; ++B) {
        BoostEdge e_ = this->base_small.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetIsomorphismMask());
      }
      for (std::tie(B, E) = large->GetEdges(); B != E; ++B) {
        BoostEdge e_ = this->base_large.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetIsomorphismMask());
      }
    }
    
    bool operator()(BoostVertex vsmall, BoostVertex vlarge) {
      return vmasks.at(vsmall) == vmasks.at(vlarge);
    }
    
    bool operator()(BoostEdge esmall, BoostEdge elarge) {
      return emasks.at(esmall) == emasks.at(elarge);
    }
    
    using BaseType::operator();
  };
  
  template <class UserCallback>
  struct MGCallback
  : public InternalCallback<IXMolecularGraph, UserCallback> {
    using BaseType = InternalCallback<IXMolecularGraph, UserCallback>;
    using Graph = typename BaseType::Graph;
    using BaseGraph = typename BaseType::BaseGraph;
    using BoostGraph = typename BaseType::BoostGraph;
    using BoostVertex = typename BaseType::BoostVertex;
    using BoostEdge = typename BaseType::BoostEdge;
    using BoostVertexIdxMap = typename BaseType::BoostVertexIdxMap;
    
    using VertexIso = Element;
    using EdgeIso = BondOrder;
    using VertexMasks = eastl::vector_map<BoostVertex, VertexIso>;
    using EdgeMasks = eastl::vector_map<BoostEdge, EdgeIso>;
    
    VertexMasks vmasks;
    EdgeMasks emasks;
    
    MGCallback(Graph small, Graph large, BoostVertexIdxMap& idxmap,
               UserCallback& callback)
    : BaseType(small, large, idxmap, callback) {
      graph::IXMolecularGraph::VertIter b, e;
      for (std::tie(b, e) = this->G_small->GetVertices(); b != e; ++b) {
        BoostVertex v = this->base_small.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetAtom()->GetElement());
      }
      for (std::tie(b, e) = this->G_large->GetVertices(); b != e; ++b) {
        BoostVertex v = this->base_large.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetAtom()->GetElement());
      }
      graph::IXMolecularGraph::EdgeIter B, E;
      for (std::tie(B, E) = this->G_small->GetEdges(); B != E; ++B) {
        BoostEdge e_ = this->base_small.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetBond()->GetOrder());
      }
      for (std::tie(B, E) = this->G_large->GetEdges(); B != E; ++B) {
        BoostEdge e_ = this->base_large.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetBond()->GetOrder());
      }
      
    }
    
    bool operator()(BoostVertex vsmall, BoostVertex vlarge) {
      return vmasks.at(vsmall) == vmasks.at(vlarge);
    }
    
    bool operator()(BoostEdge esmall, BoostEdge elarge) {
      return emasks.at(esmall) == emasks.at(elarge);
    }
    
    using BaseType::operator();
  };
  
  template <class GraphType, class InternalCallbackType, class UserCallbackType>
  void __sg_iso_runner(std::shared_ptr<GraphType> G1,
                       std::shared_ptr<GraphType> G2,
                       UserCallbackType& user_callback) {
    using namespace boost;
    using BaseGraph = typename GraphType::graph_type;
    using BoostGraph = typename BaseGraph::graph_type;
    using BoostVertex = typename BoostGraph::vertex_descriptor;
    using BoostVertexMap = eastl::vector_map<BoostVertex, size_t>;
    
    if (G2->NumVertices() < G1->NumVertices()) G1.swap(G2);
    
    // Get the graphs out
    BaseGraph& base_small = graph::access::member_graph(G1);
    BaseGraph& base_large = graph::access::member_graph(G2);
    BoostGraph& boost_small = graph::access::member_graph(base_small);
    BoostGraph& boost_large = graph::access::member_graph(base_large);
    
    // Make the vertex map
    BoostVertexMap vmap;
    size_ i = 0;
    for (auto vs = vertices(boost_small); vs.first != vs.second; ++vs.first)
      vmap.emplace(*vs.first, i++);
    i = 0;
    for (auto vs = vertices(boost_large); vs.first != vs.second; ++vs.first)
      vmap.emplace(*vs.first, i++);
    
    // Determine the order of vertices
    std::vector<BoostVertex> v_order(vertices(boost_small).first,
                                     vertices(boost_small).second);
    std::sort(v_order.begin(), v_order.end(),
              [&boost_small](BoostVertex v1, BoostVertex v2) {
                return degree(v1, boost_small) < degree(v2, boost_small);
              });
    
    // Make the callback thing
    InternalCallbackType callback(G1, G2, vmap, user_callback);
    associative_property_map<BoostVertexMap> propmap(vmap);
    
    // Run the isomorphism checking
    vf2_subgraph_iso(boost_small,  // small graph
                     boost_large,  // large graph
                     callback,     // callback object
                     propmap, // small index map
                     propmap, // large index map
                     v_order, // small vertex order,
                     callback, // edge equivalence
                     callback); // vertex equivalence
  }
  
  void SubgraphIsomorphisms(CondensedMolecularGraph G1,
                            CondensedMolecularGraph G2,
                            MappingCallback<IXCondensedMolecularGraph>& callback) {
    using UserCallback = MappingCallback<IXCondensedMolecularGraph>;
    using Graph = IXCondensedMolecularGraph;
    using Internal = CMGCallback<UserCallback>;
    __sg_iso_runner<Graph, Internal, UserCallback>(G1, G2, callback);
  }
  
  void SubgraphIsomorphisms(graph::MolecularGraph G1,
                            graph::MolecularGraph G2,
                            MappingCallback<IXMolecularGraph>& callback) {
    using UserCallback = MappingCallback<IXMolecularGraph>;
    using Graph = IXMolecularGraph;
    using Internal = MGCallback<UserCallback>;
    __sg_iso_runner<Graph, Internal, UserCallback>(G1, G2, callback);
  }
}
