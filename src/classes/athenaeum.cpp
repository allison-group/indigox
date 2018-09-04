#include <algorithm>
#include <numeric>

#include <EASTL/iterator.h>
#include <EASTL/vector_set.h>

#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/classes/athenaeum.hpp>

namespace indigox {
  using MolFrags = IXAthenaeum::MolFrags;
  // Default states for Athenaeum settings
  uint_ IXAthenaeum::Settings::AutomaticVertexLimit = 40;
  bool IXAthenaeum::Settings::AutomaticFragmentCycles = false;
  uint_ IXAthenaeum::Settings::AutomaticMaximumCycleSize = 8;
  
  IXAthenaeum::IXAthenaeum() : _overlap(0), _roverlap(0) { }
  
  IXAthenaeum::IXAthenaeum(uint_ overlap, uint_ ring_overlap)
  : _overlap(overlap), _roverlap(ring_overlap) { }
  
  bool __allow_bond_break(const Bond&) {
    return true;
  }
  
  size_ IXAthenaeum::AddAllFragments(const Molecule &mol) {
    using namespace indigox::graph;
    // Check vertex limit (throw?)
    if (mol->NumAtoms() > Settings::AutomaticVertexLimit) return 0;
    
    // Get the condensed graph
    CondensedMolecularGraph G;
    if (_graphs.find(mol) != _graphs.end()) G = _graphs.at(mol);
    else {
      G = CondenseMolecularGraph(mol->GetGraph());
      _graphs.emplace(mol, G);
      _frags.emplace(G, MolFrags());
    }
    
    eastl::vector_set<CMGVertex> all_v(G->GetVertices().first,
                                       G->GetVertices().second);
    eastl::vector_set<CMGEdge> all_e(G->GetEdges().first,
                                     G->GetEdges().second);
    
    // Get all the subgraphs
    std::vector<CondensedMolecularGraph> sub_graphs;
    algorithm::ConnectedSubgraphs(G, sub_graphs);
    size_ num = 0;
    
    // Iterate over all subgraphs
    for (CondensedMolecularGraph frag_g : sub_graphs) {
      // All vertices and edges in frag_g
      eastl::vector_set<CMGVertex> frag_v(frag_g->GetVertices().first,
                                          frag_g->GetVertices().second);
      eastl::vector_set<CMGEdge> frag_e(frag_g->GetEdges().first,
                                        frag_g->GetEdges().second);
      
      // All vertices and edges not in frag_g
      eastl::vector_set<CMGVertex> other_v;
      eastl::vector_set<CMGEdge> other_e;
      std::set_difference(all_v.begin(), all_v.end(),
                          frag_v.begin(), frag_v.end(),
                          eastl::inserter(other_v, other_v.begin()));
      std::set_difference(all_e.begin(), all_e.end(),
                          frag_e.begin(), frag_e.end(),
                          eastl::inserter(other_e, other_e.begin()));

      // Discard graph from all vertices not in f_g
      CondensedMolecularGraph d_g = G->Subgraph(other_v.begin(), other_v.end(),
                                                other_e.begin(), other_e.end());
      
      // Find the edges that have been cut
      std::vector<std::pair<CMGEdge, std::pair<CMGVertex, CMGVertex>>> cut_e;
      for (CMGEdge e : all_e) {
        if (frag_g->HasEdge(e) || d_g->HasEdge(e)) continue;
        if (!__allow_bond_break(e->GetSource()->GetBond())) {
          cut_e.clear();
          break;
        }
        CMGVertex source, target;
        std::tie(source, target) = G->GetVertices(e);
        if (frag_g->HasVertex(source)) cut_e.push_back({e, {source, target}});
        else cut_e.push_back({e, {target, source}});
      }
      if (!cut_e.size()) continue;
      
      // Find the overlap vertices from each cut
      eastl::vector_set<CMGVertex> overlaps;
      for (auto e : cut_e) {
        other_v.insert(e.second.first);
        d_g = G->Subgraph(other_v.begin(), other_v.end(),
                          other_e.begin(), other_e.end());
        for (auto vs = d_g->GetVertices(); vs.first != vs.second; ++vs.first) {
          CMGVertex v = *vs.first;
          if (v == e.second.first) continue;
          std::vector<algorithm::Path<CMGVertex>> paths;
          algorithm::AllSimplePaths(d_g, v, e.second.first, paths, _overlap);
          for (auto& path : paths) overlaps.insert(path.begin(), path.end());
        }
        other_v.erase(e.second.first);
      }
      
      // Create the fragment
      std::vector<CMGVertex> f(frag_v.begin(), frag_v.end());
      std::vector<CMGVertex> o(overlaps.begin(), overlaps.end());
      Fragment fragment = std::make_shared<IXFragment>(G, f, o);
      
      // Add the fragment if it isn't added already
      MolFrags& fs = _frags.at(G);
      if (std::find(fs.begin(), fs.end(), fragment) == fs.end()) {
        fs.emplace_back(fragment);
        ++num;
      }
    }
    return num;
  }
  
  bool IXAthenaeum::AddFragment(const Molecule &mol, Fragment frag) {
    // Get the condensed graph
    graph::CondensedMolecularGraph G;
    if (_graphs.find(mol) != _graphs.end()) G = _graphs.at(mol);
    else {
      G = graph::CondenseMolecularGraph(mol->GetGraph());
      _graphs.emplace(mol, G);
      _frags.emplace(G, MolFrags());
    }
    
    // Add the fragment
    MolFrags& fs = _frags.at(G);
    if (std::find(fs.begin(), fs.end(), frag) == fs.end()) {
      fs.emplace_back(frag);
      return true;
    }
    return false;
  }
  
  size_ IXAthenaeum::NumFragments() const {
    return std::accumulate(_frags.begin(), _frags.end(), 0,
                           [](size_ i, auto& j){ return i + j.second.size(); });
  }
  
  size_ IXAthenaeum::NumFragments(const Molecule &mol) const {
    if (_graphs.find(mol) == _graphs.end()) return 0;
    return _frags.at(_graphs.at(mol)).size();
  }
  
  const MolFrags& IXAthenaeum::GetFragments(const Molecule &mol) const {
    return _frags.at(_graphs.at(mol));
  }
  
}
