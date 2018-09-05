#include <boost/graph/vf2_sub_graph_iso.hpp>

#include <EASTL/vector_map.h>

#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>

namespace indigox::algorithm {
  
  struct CMGCallback {
    using Graph = graph::CondensedMolecularGraph;
    using BaseGraph = graph::IXCondensedMolecularGraph::graph_type;
    using BoostGraph = BaseGraph::graph_type;
    using BoostVertex = BoostGraph::vertex_descriptor;
    using BoostEdge = BoostGraph::edge_descriptor;
    using BoostVertexMap = eastl::vector_map<BoostVertex, size_>;
    using Mappings = std::map<graph::CMGVertex, std::vector<graph::CMGVertex>>;
    
    using VertexIso = graph::VertexIsoMask;
    using EdgeIso = graph::EdgeIsoMask;
    using VertexMasks = eastl::vector_map<BoostVertex, VertexIso>;
    using EdgeMasks = eastl::vector_map<BoostEdge, EdgeIso>;
    
    Graph G_small, G_large;
    BaseGraph &base_small, &base_large;
    BoostVertexMap &vmap;
    Mappings& maps;
    VertexMasks vmasks;
    EdgeMasks emasks;
    
    CMGCallback(Graph small, Graph large, BoostVertexMap& vmap2, Mappings& m)
    : G_small(small), G_large(large), base_small(graph::access::member_graph(small)),
    base_large(graph::access::member_graph(large)), vmap(vmap2), maps(m) {
      maps.clear();
      graph::IXCondensedMolecularGraph::VertIter b, e;
      for (std::tie(b, e) = G_small->GetVertices(); b != e; ++b) {
        maps.emplace(*b, std::vector<graph::CMGVertex>());
        BoostVertex v = base_small.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetIsomorphismMask());
      }
      for (std::tie(b, e) = G_large->GetVertices(); b != e; ++b) {
        BoostVertex v = base_large.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetIsomorphismMask());
      }
      graph::IXCondensedMolecularGraph::EdgeIter B, E;
      for (std::tie(B, E) = G_small->GetEdges(); B != E; ++B) {
        BoostEdge e_ = base_small.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetIsomorphismMask());
      }
      for (std::tie(B, E) = G_large->GetEdges(); B != E; ++B) {
        BoostEdge e_ = base_large.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetIsomorphismMask());
      }
      
    }
    
    template <class CMap1to2, class CMap2to1>
    bool operator()(CMap1to2 one, CMap2to1) {
      auto verts = boost::vertices(graph::access::member_graph(base_small));
      for (; verts.first != verts.second; ++verts.first) {
        auto vsmall = base_small.GetVertex(*verts.first);
        auto vlarge = base_large.GetVertex(boost::get(one, *verts.first));
        maps.at(vsmall->shared_from_this()).push_back(vlarge->shared_from_this());
      }
      return true;
    }
    
    bool operator()(BoostVertex vsmall, BoostVertex vlarge) {
      return vmasks.at(vsmall) == vmasks.at(vlarge);
    }
    
    bool operator()(BoostEdge esmall, BoostEdge elarge) {
      return emasks.at(esmall) == emasks.at(elarge);
    }
  };
  
  struct MGCallback {
    using Graph = graph::MolecularGraph;
    using BaseGraph = graph::IXMolecularGraph::graph_type;
    using BoostGraph = BaseGraph::graph_type;
    using BoostVertex = BoostGraph::vertex_descriptor;
    using BoostEdge = BoostGraph::edge_descriptor;
    using BoostVertexMap = eastl::vector_map<BoostVertex, size_>;
    using Mappings = std::map<graph::MGVertex, std::vector<graph::MGVertex>>;
    
    using VertexIso = Element;
    using EdgeIso = BondOrder;
    using VertexMasks = eastl::vector_map<BoostVertex, VertexIso>;
    using EdgeMasks = eastl::vector_map<BoostEdge, EdgeIso>;
    
    Graph G_small, G_large;
    BaseGraph &base_small, &base_large;
    BoostVertexMap &vmap;
    Mappings& maps;
    VertexMasks vmasks;
    EdgeMasks emasks;
    
    MGCallback(Graph small, Graph large, BoostVertexMap& vmap2, Mappings& m)
    : G_small(small), G_large(large), base_small(graph::access::member_graph(small)),
    base_large(graph::access::member_graph(large)), vmap(vmap2), maps(m) {
      maps.clear();
      graph::IXMolecularGraph::VertIter b, e;
      for (std::tie(b, e) = G_small->GetVertices(); b != e; ++b) {
        maps.emplace(*b, std::vector<graph::MGVertex>());
        BoostVertex v = base_small.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetAtom()->GetElement());
      }
      for (std::tie(b, e) = G_large->GetVertices(); b != e; ++b) {
        BoostVertex v = base_large.GetDescriptor(b->get());
        vmasks.emplace(v, (*b)->GetAtom()->GetElement());
      }
      graph::IXMolecularGraph::EdgeIter B, E;
      for (std::tie(B, E) = G_small->GetEdges(); B != E; ++B) {
        BoostEdge e_ = base_small.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetBond()->GetOrder());
      }
      for (std::tie(B, E) = G_large->GetEdges(); B != E; ++B) {
        BoostEdge e_ = base_large.GetDescriptor(B->get());
        emasks.emplace(e_, (*B)->GetBond()->GetOrder());
      }
      
    }
    
    template <class CMap1to2, class CMap2to1>
    bool operator()(CMap1to2 one, CMap2to1) {
      auto verts = boost::vertices(graph::access::member_graph(base_small));
      for (; verts.first != verts.second; ++verts.first) {
        auto vsmall = base_small.GetVertex(*verts.first);
        auto vlarge = base_large.GetVertex(boost::get(one, *verts.first));
        maps.at(vsmall->shared_from_this()).push_back(vlarge->shared_from_this());
      }
      return true;
    }
    
    bool operator()(BoostVertex vsmall, BoostVertex vlarge) {
      return vmasks.at(vsmall) == vmasks.at(vlarge);
    }
    
    bool operator()(BoostEdge esmall, BoostEdge elarge) {
      return emasks.at(esmall) == emasks.at(elarge);
    }
  };
  
  template <class GraphType, class CallbackType>
  void __sg_iso_runner(std::shared_ptr<GraphType> G1,
                       std::shared_ptr<GraphType> G2,
                       std::map<typename GraphType::VertexType,
                             std::vector<typename GraphType::VertexType>>& maps) {
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
    CallbackType callback(G1, G2, vmap, maps);
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
  
  void SubgraphIsomorphisms(graph::CondensedMolecularGraph G1,
                            graph::CondensedMolecularGraph G2,
                            std::map<graph::CMGVertex,
                            std::vector<graph::CMGVertex>>& maps) {
    __sg_iso_runner<graph::IXCondensedMolecularGraph, CMGCallback>(G1, G2, maps);
  }
  
  void SubgraphIsomorphisms(graph::MolecularGraph G1,
                            graph::MolecularGraph G2,
                            std::map<graph::MGVertex,
                            std::vector<graph::MGVertex>>& maps) {
    __sg_iso_runner<graph::IXMolecularGraph, MGCallback>(G1, G2, maps);
  }
}
