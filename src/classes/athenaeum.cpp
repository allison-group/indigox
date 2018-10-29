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
  // ===========================================================================
  // == Fragment implementation ================================================
  // ===========================================================================
  
  struct Fragment::FragmentData {
    graph::sCondensedMolecularGraph graph;
    std::vector<graph::CMGVertex> frag;
    std::vector<Fragment::OverlapVertex> overlap;
    std::vector<Fragment::AtmType> atoms;
    std::vector<Fragment::BndType> bonds;
    std::vector<Fragment::AngType> angles;
    std::vector<Fragment::DhdType> dihedrals;
  };
  
  void _MGVertexToCMGVertex(std::vector<graph::MGVertex>& v_in,
                            std::vector<graph::CMGVertex>& v_out,
                            graph::CondensedMolecularGraph& g) {
    v_out.clear();
    eastl::vector_set<graph::MGVertex> contracted, relaxed;
    for (graph::MGVertex &v : v_in) {
      if (g.HasVertex(v)) relaxed.emplace(v);
      else if (g.HasCondensedVertex(v)) contracted.emplace(v);
      else throw std::runtime_error("Unknown vertex provided");
    }
    for (graph::MGVertex v : relaxed) {
      graph::CMGVertex v_con = g.GetVertex(v);
      for (graph::MGVertex &u : v_con.GetContractedVertices()) {
        if (contracted.find(u) == contracted.end())
          throw std::runtime_error("Missing a contracted vertex from input");
        contracted.erase(u);
      }
      v_out.emplace_back(v_con);
    }
    if (!contracted.empty())
      throw std::runtime_error("Not all required vertices provided");
  }
  
  Fragment::OverlapType _GetOverType(const graph::CMGVertex& v,
                                     graph::CondensedMolecularGraph& g) {
    if (!g.HasVertex(v)) throw std::runtime_error("Issue with overlap type");
    return Fragment::OverlapType::GenericOverlap;
  }
  
  Fragment::Fragment() : _dat(nullptr) { }
  Fragment::Fragment(const Fragment& frag) : _dat(frag._dat) { }
  Fragment::Fragment(Fragment&& frag) : _dat(std::move(frag._dat)) { }
  Fragment& Fragment::operator=(const Fragment &frag) {
    if (&frag != this) _dat = frag._dat;
    return *this;
  }
  Fragment& Fragment::operator=(Fragment &&frag) {
    _dat = std::move(frag._dat);
    return *this;
  }
  
  Fragment::Fragment(graph::MolecularGraph& G,
                     std::vector<graph::MGVertex>& frag,
                     std::vector<graph::MGVertex>& overlap)
  : _dat(std::make_shared<FragmentData>()) {
    if (frag.empty()) throw std::runtime_error("A fragment needs vertices");
    // Induce a new subgraph
    graph::CondensedMolecularGraph& CG = G.GetCondensedGraph();
    std::vector<graph::CMGVertex> contracted_overlap, combined;
    _MGVertexToCMGVertex(frag, _dat->frag, CG);
    _MGVertexToCMGVertex(overlap, contracted_overlap, CG);
    combined.assign(_dat->frag.begin(),
                    _dat->frag.end());
    combined.insert(combined.end(), contracted_overlap.begin(),
                    contracted_overlap.end());
    _dat->graph = CG.Subgraph(combined);
    _dat->atoms.assign(frag.begin(), frag.end());
    _dat->overlap.reserve(contracted_overlap.size());
    for (graph::CMGVertex& v : contracted_overlap)
      _dat->overlap.emplace_back(_GetOverType(v, *_dat->graph), v);
    
    // Combine all MGVertices for checking purposes
    std::vector<AtmType> atm_check(frag.begin(), frag.end());
    atm_check.insert(atm_check.end(), overlap.begin(), overlap.end());
    size_t acceptable_pos = frag.size();
    
    // Get the molecule
    Molecule& mol = frag.front().GetAtom().GetMolecule();
    
    // Get the bonds which are allowed
    std::deque<BndType> tmp_bnd;
    for (sBond bnd : mol.GetBonds()) {
      AtmType a = G.GetVertex(bnd->GetSourceAtom());
      AtmType b = G.GetVertex(bnd->GetTargetAtom());
      auto a_pos = std::find(atm_check.begin(), atm_check.end(), a);
      if (a_pos == atm_check.end()) continue;
      auto b_pos = std::find(atm_check.begin(), atm_check.end(), b);
      if (b_pos == atm_check.end()) continue;
      size_t dangle = 0;
      if (std::distance(atm_check.begin(), a_pos) >= acceptable_pos) ++dangle;
      if (std::distance(atm_check.begin(), b_pos) >= acceptable_pos) ++dangle;
      // Bonds are allowed one dangling. Danglings go on end
      if (dangle > 1) continue;
      if (dangle) tmp_bnd.emplace_back(a, b);
      else tmp_bnd.emplace_front(a, b);
    }
    _dat->bonds.assign(tmp_bnd.begin(), tmp_bnd.end());
    
    // Get the angles which are allowed
    std::deque<AngType> tmp_ang;
    for (sAngle ang : mol.GetAngles()) {
      AtmType a = G.GetVertex(ang->GetAtoms().first);
      AtmType b = G.GetVertex(ang->GetAtoms().second);
      AtmType c = G.GetVertex(ang->GetAtoms().third);
      auto a_pos = std::find(atm_check.begin(), atm_check.end(), a);
      if (a_pos == atm_check.end()) continue;
      auto b_pos = std::find(atm_check.begin(), atm_check.end(), b);
      if (b_pos == atm_check.end()) continue;
      auto c_pos = std::find(atm_check.begin(), atm_check.end(), c);
      if (c_pos == atm_check.end()) continue;
      size_t dangle = 0;
      if (std::distance(atm_check.begin(), a_pos) >= acceptable_pos) ++dangle;
      if (std::distance(atm_check.begin(), b_pos) >= acceptable_pos) ++dangle;
      if (std::distance(atm_check.begin(), c_pos) >= acceptable_pos) ++dangle;
      // Angles are allowed one dangling. Danglings go on end
      if (dangle > 1) continue;
      if (dangle) tmp_ang.emplace_back(a, b, c);
      else tmp_ang.emplace_front(a, b, c);
    }
    _dat->angles.assign(tmp_ang.begin(), tmp_ang.end());
    
    // Get the dihedrals which are allowed
    std::deque<DhdType> tmp_dhd;
    for (sDihedral dhd : mol.GetDihedrals()) {
      AtmType a = G.GetVertex(dhd->GetAtoms().first);
      AtmType b = G.GetVertex(dhd->GetAtoms().second);
      AtmType c = G.GetVertex(dhd->GetAtoms().third);
      AtmType d = G.GetVertex(dhd->GetAtoms().fourth);
      auto a_pos = std::find(atm_check.begin(), atm_check.end(), a);
      if (a_pos == atm_check.end()) continue;
      auto b_pos = std::find(atm_check.begin(), atm_check.end(), b);
      if (b_pos == atm_check.end()) continue;
      auto c_pos = std::find(atm_check.begin(), atm_check.end(), c);
      if (c_pos == atm_check.end()) continue;
      auto d_pos = std::find(atm_check.begin(), atm_check.end(), d);
      if (d_pos == atm_check.end()) continue;
      bool d1 = std::distance(atms_check.begin(), a_pos) >= acceptable_pos;
      bool d2 = std::distance(atms_check.begin(), b_pos) >= acceptable_pos;
      bool d3 = std::distance(atms_check.begin(), c_pos) >= acceptable_pos;
      bool d4 = std::distance(atms_check.begin(), d_pos) >= acceptable_pos;
      // Need two adjacent atoms to not be dangling
      if (!((d1 && d2) || (d2 && d3) || (d3 && d4))) continue;
      if (d1 || d2 || d3 || d4) tmp_dhd.emplace_back(a, b, c, d);
      else tmp_dhd.emplace_front(a, b, c, d);
    }
    _dat->dihedrals.assign(tmp_dhd.begin(), tmp_dhd.end());
  }
  
  graph::CondensedMolecularGraph& Fragment::GetGraph() const {
    return *_dat->graph;
  }
  
  const std::vector<graph::CMGVertex>& Fragment::GetFragment() const {
    return _dat->frag;
  }
  
  const std::vector<Fragment::OverlapVertex>& Fragment::GetOverlap() const {
    return _dat->overlap;
  }
  
  bool Fragment::IsFragmentVertex(graph::CMGVertex &v) const {
    return (std::find(_dat->frag.begin(), _dat->frag.end(), v)
            != _dat->frag.end());
  }
  
  bool Fragment::IsOverlapVertex(graph::CMGVertex &v) const {
    return (std::find(_dat->overlap.begin(), _dat->overlap.end(), v)
            != _dat->overlap.end());
  }
  
  const std::vector<Fragment::AtmType> Fragment::GetAtoms() const {
    return _dat->atoms;
  }
  
  const std::vector<Fragment::BndType> Fragment::GetBonds() const {
    return _dat->bonds;
  }
  
  const std::vector<Fragment::AngType> Fragment::GetAngles() const {
    return _dat->angles;
  }
  
  const std::vector<Fragment::DhdType> Fragment::GetDihedrals() const {
    return _dat->dihedrals;
  }
  
}

//namespace indigox {
//  using MolFrags = IXAthenaeum::MolFrags;
//  // Default states for Athenaeum settings
//  uint32_t IXAthenaeum::Settings::AutomaticVertexLimit = 40;
//  bool IXAthenaeum::Settings::AutomaticFragmentCycles = false;
//  uint32_t IXAthenaeum::Settings::AutomaticMaximumCycleSize = 8;
//  
//  
//  
//  IXAthenaeum::IXAthenaeum(Forcefield ff) : indigox::IXAthenaeum(ff, 0, 0) { }
//  
//  IXAthenaeum::IXAthenaeum(Forcefield ff, uint32_t overlap, uint32_t ring_overlap)
//  : _ff(ff), _overlap(overlap), _roverlap(ring_overlap), _man(false) { }
//  
//  bool __allow_bond_break(const Bond&) {
//    return true;
//  }
//  
//  size_t IXAthenaeum::AddAllFragments(const Molecule &mol) {
//    using namespace indigox::graph;
//    // Check vertex limit (throw?)
//    if (mol->NumAtoms() > Settings::AutomaticVertexLimit) return 0;
//    
//    // Get the condensed graph
//    CondensedMolecularGraph G;
//    if (_graphs.find(mol) != _graphs.end()) G = _graphs.at(mol);
//    else {
//      G = CondenseMolecularGraph(mol->GetGraph());
//      _graphs.emplace(mol, G);
//      _frags.emplace(G, MolFrags());
//    }
//    
//    eastl::vector_set<CMGVertex> all_v(G->GetVertices().first,
//                                       G->GetVertices().second);
//    eastl::vector_set<CMGEdge> all_e(G->GetEdges().first,
//                                     G->GetEdges().second);
//    
//    // Get all the subgraphs
//    std::vector<CondensedMolecularGraph> sub_graphs;
//    algorithm::ConnectedSubgraphs(G, sub_graphs);
//    size_t num = 0;
//    
//    // Iterate over all subgraphs
//    for (CondensedMolecularGraph frag_g : sub_graphs) {
//      // All vertices and edges in frag_g
//      eastl::vector_set<CMGVertex> frag_v(frag_g->GetVertices().first,
//                                          frag_g->GetVertices().second);
//      eastl::vector_set<CMGEdge> frag_e(frag_g->GetEdges().first,
//                                        frag_g->GetEdges().second);
//      
//      // All vertices and edges not in frag_g
//      eastl::vector_set<CMGVertex> other_v;
//      eastl::vector_set<CMGEdge> other_e;
//      std::set_difference(all_v.begin(), all_v.end(),
//                          frag_v.begin(), frag_v.end(),
//                          eastl::inserter(other_v, other_v.begin()));
//      std::set_difference(all_e.begin(), all_e.end(),
//                          frag_e.begin(), frag_e.end(),
//                          eastl::inserter(other_e, other_e.begin()));
//
//      // Discard graph from all vertices not in f_g
//      CondensedMolecularGraph d_g = G->Subgraph(other_v.begin(), other_v.end(),
//                                                other_e.begin(), other_e.end());
//      
//      // Find the edges that have been cut
//      std::vector<std::pair<CMGEdge, std::pair<CMGVertex, CMGVertex>>> cut_e;
//      for (CMGEdge e : all_e) {
//        if (frag_g->HasEdge(e) || d_g->HasEdge(e)) continue;
//        if (!__allow_bond_break(e->GetSource()->GetBond())) {
//          cut_e.clear();
//          break;
//        }
//        CMGVertex source, target;
//        std::tie(source, target) = G->GetVertices(e);
//        if (frag_g->HasVertex(source)) cut_e.push_back({e, {source, target}});
//        else cut_e.push_back({e, {target, source}});
//      }
//      if (!cut_e.size()) continue;
//      
//      // Find the overlap vertices from each cut
//      eastl::vector_set<CMGVertex> overlaps;
//      for (auto e : cut_e) {
//        other_v.insert(e.second.first);
//        d_g = G->Subgraph(other_v.begin(), other_v.end(),
//                          other_e.begin(), other_e.end());
//        for (auto vs = d_g->GetVertices(); vs.first != vs.second; ++vs.first) {
//          CMGVertex v = *vs.first;
//          if (v == e.second.first) continue;
//          std::vector<algorithm::Path<CMGVertex>> paths;
//          algorithm::AllSimplePaths(d_g, v, e.second.first, paths, _overlap);
//          for (auto& path : paths) overlaps.insert(path.begin(), path.end());
//        }
//        other_v.erase(e.second.first);
//      }
//      
//      // Create the fragment
//      std::vector<CMGVertex> f(frag_v.begin(), frag_v.end());
//      std::vector<CMGVertex> o(overlaps.begin(), overlaps.end());
//      Fragment fragment = std::make_shared<IXFragment>(G, mol, f, o);
//      
//      // Add the fragment if it isn't added already
//      MolFrags& fs = _frags.at(G);
//      if (std::find(fs.begin(), fs.end(), fragment) == fs.end()) {
//        fs.emplace_back(fragment);
//        ++num;
//      }
//    }
//    return num;
//  }
//  
//  bool IXAthenaeum::AddFragment(const Molecule &mol, Fragment frag) {
//    // Get the condensed graph
//    graph::CondensedMolecularGraph G;
//    if (_graphs.find(mol) != _graphs.end()) G = _graphs.at(mol);
//    else {
//      G = graph::CondenseMolecularGraph(mol->GetGraph());
//      _graphs.emplace(mol, G);
//      _frags.emplace(G, MolFrags());
//    }
//    
//    // Add the fragment
//    MolFrags& fs = _frags.at(G);
//    if (std::find(fs.begin(), fs.end(), frag) == fs.end()) {
//      fs.emplace_back(frag);
//      return true;
//    }
//    return false;
//  }
//  
//  size_t IXAthenaeum::NumFragments() const {
//    return std::accumulate(_frags.begin(), _frags.end(), 0,
//                           [](size_t i, auto& j){ return i + j.second.size(); });
//  }
//  
//  size_t IXAthenaeum::NumFragments(const Molecule &mol) const {
//    if (_graphs.find(mol) == _graphs.end()) return 0;
//    return _frags.at(_graphs.at(mol)).size();
//  }
//  
//  const MolFrags& IXAthenaeum::GetFragments(const Molecule &mol) const {
//    return _frags.at(_graphs.at(mol));
//  }
//  
//}

