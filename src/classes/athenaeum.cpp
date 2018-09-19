#include <algorithm>
#include <deque>
#include <numeric>

#include <EASTL/iterator.h>
#include <EASTL/vector_set.h>

#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/angle.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>


namespace indigox {
  using MolFrags = IXAthenaeum::MolFrags;
  // Default states for Athenaeum settings
  uint_ IXAthenaeum::Settings::AutomaticVertexLimit = 40;
  bool IXAthenaeum::Settings::AutomaticFragmentCycles = false;
  uint_ IXAthenaeum::Settings::AutomaticMaximumCycleSize = 8;
  
  IXFragment::IXFragment(const graph::CondensedMolecularGraph &G,
                         const Molecule &mol,
                         const std::vector<graph::CMGVertex> &frag,
                         const std::vector<graph::CMGVertex> &overlap)
  : _mol(mol), _frag(frag.begin(), frag.end()),
  _overlap(overlap.begin(), overlap.end()) {
    std::vector<graph::CMGVertex> tmp(frag.begin(), frag.end());
    tmp.insert(tmp.end(), overlap.begin(), overlap.end());
    _g = G->InduceSubgraph(tmp.begin(), tmp.end());
    graph::MolecularGraph molG = _g->GetSource();
    
    std::vector<Vert> atms_check;
    // Get the MGVertices that are part of the fragment/overlap
    for (graph::CMGVertex v : frag) {
      _atms.emplace_back(v->GetSource());
      for (auto& cv : v->GetContractedVertices()) _atms.emplace_back(cv.second);
    }
    int d_max = std::distance(_atms.begin(), _atms.end());
    atms_check.assign(_atms.begin(), _atms.end());
    
    for (graph::CMGVertex v : overlap) {
      atms_check.emplace_back(v->GetSource());
      for (auto& cv : v->GetContractedVertices())
        atms_check.emplace_back(cv.second);
    }
    
    std::deque<std::pair<Vert, Vert>> tmp_bnd;
    // Get the bonds which are allowed, including dangling bonds
    for (auto it = mol->GetBonds(); it.first != it.second; ++it.first) {
      auto atms = (*it.first)->GetAtoms();
      graph::MGVertex v1 = molG->GetVertex(atms.first);
      graph::MGVertex v2 = molG->GetVertex(atms.second);
      auto v1_pos = std::find(atms_check.begin(), atms_check.end(), v1);
      if (v1_pos == atms_check.end()) continue;
      auto v2_pos = std::find(atms_check.begin(), atms_check.end(), v2);
      if (v2_pos == atms_check.end()) continue;
      // Only one atom of the bond can be dangling
      size_t dangle_count = 0;
      if (std::distance(atms_check.begin(), v1_pos) >= d_max) ++dangle_count;
      if (std::distance(atms_check.begin(), v2_pos) >= d_max) ++dangle_count;
      if (dangle_count > 1) continue;
      if (dangle_count) tmp_bnd.emplace_back(v1,v2);
      else tmp_bnd.emplace_front(v1,v2);
    }
    _bnds.assign(tmp_bnd.begin(), tmp_bnd.end());
    
    std::deque<stdx::triple<Vert, Vert, Vert>> tmp_ang;
    for (auto it = mol->GetAngles(); it.first != it.second; ++it.first) {
      auto atms = (*it.first)->GetAtoms();
      graph::MGVertex v1 = molG->GetVertex(atms.first);
      graph::MGVertex v2 = molG->GetVertex(atms.second);
      graph::MGVertex v3 = molG->GetVertex(atms.third);
      auto v1_pos = std::find(atms_check.begin(), atms_check.end(), v1);
      if (v1_pos == atms_check.end()) continue;
      auto v2_pos = std::find(atms_check.begin(), atms_check.end(), v2);
      if (v2_pos == atms_check.end()) continue;
      auto v3_pos = std::find(atms_check.begin(), atms_check.end(), v3);
      if (v3_pos == atms_check.end()) continue;
      // Only one atom of the angle can be dangling
      size_t dangle_count = 0;
      if (std::distance(atms_check.begin(), v1_pos) >= d_max) ++dangle_count;
      if (std::distance(atms_check.begin(), v2_pos) >= d_max) ++dangle_count;
      if (std::distance(atms_check.begin(), v3_pos) >= d_max) ++dangle_count;
      if (dangle_count > 1) continue;
      if (dangle_count) tmp_ang.emplace_back(v1,v2,v3);
      else tmp_ang.emplace_front(v1,v2,v3);
    }
    _angs.assign(tmp_ang.begin(), tmp_ang.end());
    
    std::deque<stdx::quad<Vert, Vert, Vert, Vert>> tmp_dhd;
    for (auto it = mol->GetDihedrals(); it.first != it.second; ++it.first) {
      auto atms = (*it.first)->GetAtoms();
      graph::MGVertex v1 = molG->GetVertex(atms.first);
      graph::MGVertex v2 = molG->GetVertex(atms.second);
      graph::MGVertex v3 = molG->GetVertex(atms.third);
      graph::MGVertex v4 = molG->GetVertex(atms.fourth);
      auto v1_pos = std::find(atms_check.begin(), atms_check.end(), v1);
      if (v1_pos == atms_check.end()) continue;
      auto v2_pos = std::find(atms_check.begin(), atms_check.end(), v2);
      if (v2_pos == atms_check.end()) continue;
      auto v3_pos = std::find(atms_check.begin(), atms_check.end(), v3);
      if (v3_pos == atms_check.end()) continue;
      auto v4_pos = std::find(atms_check.begin(), atms_check.end(), v4);
      if (v4_pos == atms_check.end()) continue;
      
      bool dangling = false;
      int d1 = std::distance(atms_check.begin(), v1_pos);
      int d2 = std::distance(atms_check.begin(), v2_pos);
      int d3 = std::distance(atms_check.begin(), v3_pos);
      int d4 = std::distance(atms_check.begin(), v4_pos);
      if (d1 >= d_max || d2 >= d_max || d3 >= d_max || d4 >= d_max)
        dangling = true;
      // Need two adjacent atoms to not be dangling
      if (!((d1 < d_max && d2 < d_max) || (d2 < d_max && d3 < d_max)
            || (d3 < d_max && d4 < d_max))) continue;
      if (dangling) tmp_dhd.emplace_back(v1,v2,v3,v4);
      else tmp_dhd.emplace_front(v1,v2,v3,v4);
    }
    _dhds.assign(tmp_dhd.begin(), tmp_dhd.end());
    
  }
  
  IXAthenaeum::IXAthenaeum(Forcefield ff) : indigox::IXAthenaeum(ff, 0, 0) { }
  
  IXAthenaeum::IXAthenaeum(Forcefield ff, uint_ overlap, uint_ ring_overlap)
  : _ff(ff), _overlap(overlap), _roverlap(ring_overlap), _man(false) { }
  
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
      Fragment fragment = std::make_shared<IXFragment>(G, mol, f, o);
      
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

