#include <algorithm>
#include <deque>
#include <iterator>
#include <numeric>
#include <vector>

#include <EASTL/iterator.h>
#include <EASTL/vector_set.h>

#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/angle.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
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
                                     graph::CondensedMolecularGraph& g,
                                     std::vector<graph::CMGVertex>&) {
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
    std::sort(_dat->frag.begin(), _dat->frag.end());
    _MGVertexToCMGVertex(overlap, contracted_overlap, CG);
    combined.assign(_dat->frag.begin(),
                    _dat->frag.end());
    combined.insert(combined.end(), contracted_overlap.begin(),
                    contracted_overlap.end());
    _dat->graph = CG.Subgraph(combined);
    if (!_dat->graph->IsConnected())
      throw std::runtime_error("A fragment must be connected");
    _dat->atoms.assign(frag.begin(), frag.end());
    _dat->overlap.reserve(contracted_overlap.size());
    for (graph::CMGVertex& v : contracted_overlap)
      _dat->overlap.emplace_back(_GetOverType(v, *_dat->graph, _dat->frag), v);
    std::sort(_dat->overlap.begin(), _dat->overlap.end());
    
    // Combine all MGVertices for checking purposes
    std::vector<AtmType> atm_check(frag.begin(), frag.end());
    atm_check.insert(atm_check.end(), overlap.begin(), overlap.end());
    long acceptable_pos = frag.size();
    
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
      bool d1 = std::distance(atm_check.begin(), a_pos) >= acceptable_pos;
      bool d2 = std::distance(atm_check.begin(), b_pos) >= acceptable_pos;
      bool d3 = std::distance(atm_check.begin(), c_pos) >= acceptable_pos;
      bool d4 = std::distance(atm_check.begin(), d_pos) >= acceptable_pos;
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
  
  size_t Fragment::Size() const { return _dat->frag.size(); }
  
  const std::vector<Fragment::OverlapVertex>& Fragment::GetOverlap() const {
    return _dat->overlap;
  }
  
  bool Fragment::IsFragmentVertex(graph::CMGVertex &v) const {
    return (std::find(_dat->frag.begin(), _dat->frag.end(), v)
            != _dat->frag.end());
  }
  
  bool Fragment::IsOverlapVertex(graph::CMGVertex &v) const {
    auto pos = std::find_if(_dat->overlap.begin(), _dat->overlap.end(),
                            [&v](auto& u) { return u.second == v; });
    return pos != _dat->overlap.end();
  }
  
  const std::vector<Fragment::AtmType>& Fragment::GetAtoms() const {
    return _dat->atoms;
  }
  
  const std::vector<Fragment::BndType>& Fragment::GetBonds() const {
    return _dat->bonds;
  }
  
  const std::vector<Fragment::AngType>& Fragment::GetAngles() const {
    return _dat->angles;
  }
  
  const std::vector<Fragment::DhdType>& Fragment::GetDihedrals() const {
    return _dat->dihedrals;
  }
  
  bool Fragment::operator==(const Fragment& frag) const {
    if (_dat == frag._dat) return true;
    if (_dat->frag != frag._dat->frag) return false;
    if (_dat->overlap != frag._dat->overlap) return false;
    return true;
  }
  
  bool Fragment::operator<(const Fragment &frag) const {
    if (_dat == frag._dat) return false;
    if (_dat->frag < frag._dat->frag) return true;
    return _dat->overlap < frag._dat->overlap;
  }
  
  bool Fragment::operator>(const Fragment &frag) const {
    if (_dat == frag._dat) return false;
    if (_dat->frag > frag._dat->frag) return true;
    return _dat->overlap > frag._dat->overlap;
  }
 
  // ===========================================================================
  // == Athenaeum implementation ===============================================
  // ===========================================================================
  
  // Default settings
  uint32_t Athenaeum::Settings::AtomLimit = 40;
  uint32_t Athenaeum::Settings::MinimumFragmentSize = 1;
  uint32_t Athenaeum::Settings::MaximumFragmentSize = 40;
  bool Athenaeum::Settings::FragmentCycles = false;
  uint32_t Athenaeum::Settings::MaximumCycleSize = 8;
  uint32_t Athenaeum::Settings::DefaultOverlap = 2;
  uint32_t Athenaeum::Settings::DefaultCycleOverlap = 2;
  
  Athenaeum::Athenaeum(Forcefield& ff)
  : Athenaeum(ff, Settings::DefaultOverlap, Settings::DefaultCycleOverlap) { }
  
  Athenaeum::Athenaeum(Forcefield& ff, uint32_t overlap)
  : Athenaeum(ff, overlap, Settings::DefaultCycleOverlap) { }
  
  Athenaeum::Athenaeum(Forcefield& ff, uint32_t overlap, uint32_t cycleoverlap)
  : _ff(ff), _overlap(overlap), _roverlap(cycleoverlap),
  _man(false), _frags() { }
  
  size_t Athenaeum::NumFragments() const {
    return std::accumulate(_frags.begin(), _frags.end(), 0,
                           [](size_t i, auto& j){ return i + j.second.size(); });
  }
  
  size_t Athenaeum::NumFragments(Molecule &mol) const {
    sMolecule m = mol.shared_from_this();
    if (_frags.find(m) == _frags.end()) return 0;
    return _frags.at(m).size();
  }
  
  const Athenaeum::MoleculeFragments& Athenaeum::GetFragments() const {
    return _frags;
  }
  
  const Athenaeum::FragContain& Athenaeum::GetFragments(Molecule& mol) const {
    sMolecule m = mol.shared_from_this();
    if (_frags.find(m) == _frags.end())
      throw std::runtime_error("No fragmenst for molecule available");
    return _frags.at(m);
  }
  
  bool Athenaeum::HasFragments(Molecule &mol) const {
    sMolecule m = mol.shared_from_this();
    return _frags.find(m) != _frags.end();
  }

  bool Athenaeum::AddFragment(Molecule &mol, Fragment &frag) {
    // Check that the fragment matches the molecule
    graph::MolecularGraph& MG = mol.GetGraph();
    graph::CondensedMolecularGraph& CG = MG.GetCondensedGraph();
    graph::CondensedMolecularGraph& fg = frag.GetGraph();
    while (fg.IsSubgraph()) { fg = fg.GetSuperGraph(); }
    if (&CG != &fg) return false;
    
    // Check that the molecule forcefield matchs the athenaeum forcefield
    if (!mol.HasForcefield()) return false;
    if (mol.GetForcefield() != _ff) return false;
    
    sMolecule m = mol.shared_from_this();
    if (_frags.find(m) == _frags.end()) _frags.emplace(m, FragContain());
    _frags.at(m).emplace_back(frag);
    return true;
  }
  
  bool CanCutEdge(graph::CMGEdge& e, graph::CondensedMolecularGraph& g) {
    return g.HasEdge(e);
  }
  
  size_t Athenaeum::AddAllFragments(Molecule& mol) {
    using namespace indigox::graph;
    // Perform checks
    if (!mol.IsFrozen())
      throw std::runtime_error("Can only add fragments from frozen molecule");
    if (!mol.HasForcefield())
      throw std::runtime_error("Attempting to fragment unparameterised molecule");
    if (mol.GetForcefield() != _ff)
      throw std::runtime_error("Forcefield mismatch");
    if (mol.NumAtoms() > Settings::AtomLimit)
      throw std::runtime_error("Molecule too large to automagically fragment");
    
    sMolecule m = mol.shared_from_this();
    if (_frags.find(m) == _frags.end()) _frags.emplace(m, FragContain());
    size_t initial_count = _frags.at(m).size();
    
    // Get all the subgraphs of the molecule's condensed graph
    MolecularGraph& MG = mol.GetGraph();
    CondensedMolecularGraph& CG = MG.GetCondensedGraph();
    std::vector<sCondensedMolecularGraph> sub_graphs;
    algorithm::ConnectedSubgraphs(CG, sub_graphs);
    eastl::vector_set<CMGVertex> all_vertices(CG.GetVertices().begin(),
                                              CG.GetVertices().end());
    eastl::vector_set<CMGEdge> all_edges(CG.GetEdges().begin(),
                                         CG.GetEdges().end());
    
    // Decide if each subgraph can be made into a fragment
    for (sCondensedMolecularGraph sub : sub_graphs) {
      // Sort the vertices/edges of CG into not in sub and in sub
      eastl::vector_set<CMGVertex> sub_vertices(sub->GetVertices().begin(),
                                                sub->GetVertices().end());
      eastl::vector_set<CMGEdge> sub_edges(sub->GetEdges().begin(),
                                           sub->GetEdges().end());
      std::vector<CMGVertex> other_vertices;
      std::vector<CMGEdge> other_edges;
      std::set_difference(all_vertices.begin(), all_vertices.end(),
                          sub_vertices.begin(), sub_vertices.end(),
                          std::back_inserter(other_vertices));
      std::set_difference(all_edges.begin(), all_edges.end(),
                          sub_edges.begin(), sub_edges.end(),
                          std::back_inserter(other_edges));
      
      // Determine which edges are cut
      // cut_edge.second is true if source vertex in other_vertices, false if not
      std::vector<std::pair<CMGEdge, bool>> cut_edges;
      for (CMGEdge e : other_edges) {
        CMGVertex u = CG.GetSourceVertex(e);
        CMGVertex v = CG.GetTargetVertex(e);
        bool has_u = sub->HasVertex(u);
        bool has_v = sub->HasVertex(v);
        if (has_u && has_v) throw std::runtime_error("WTF?!");
        else if (has_u && !has_v) cut_edges.emplace_back(e, false);
        else if (!has_u && has_v) cut_edges.emplace_back(e, true);
        else continue;
        // May as well check cutabliity of edge at same time
        if (!CanCutEdge(e, CG)) {
          cut_edges.clear();
          break;
        }
      }
      if (cut_edges.empty()) continue;
      
      // Find all the vertices within _overlap of the fragment vertices
      eastl::vector_set<CMGVertex> overlap_vertices;
      for (CMGVertex v : sub_vertices) {
        for (CMGVertex u : other_vertices) {
          if (overlap_vertices.find(u) != overlap_vertices.end()) continue;
          auto path = algorithm::ShortestPath(CG, u, v);
          if (!path.empty() && path.size() <= _overlap)
            overlap_vertices.emplace(u);
        }
      }
      
      // every leaf in overlap must have minimum path length of _overlap to
      // each vertex in fragment
      std::vector<CMGVertex> fragoververt(sub_vertices.begin(),
                                          sub_vertices.end());
      fragoververt.insert(fragoververt.end(), overlap_vertices.begin(),
                          overlap_vertices.end());
      sCondensedMolecularGraph withoverlap = CG.Subgraph(fragoververt);
      CondensedMolecularGraph::ComponentContain tmp;
      if (algorithm::ConnectedComponents(*withoverlap, tmp) > 1) continue;
      bool bad_overlaps = false;
      for (CMGVertex u : overlap_vertices) {
        if (withoverlap->Degree(u) > 1) continue;
        for (CMGVertex v : sub_vertices) {
          auto path = algorithm::ShortestPath(*withoverlap, u, v);
          if (path.size() < _overlap) bad_overlaps = true;
        }
      }
      if (bad_overlaps) continue;
      
      // Create the fragment and add it in
      std::vector<MGVertex> final_frag, final_overlap;
      for (CMGVertex v : sub_vertices) {
        final_frag.emplace_back(v.GetSource());
        auto contract = v.GetContractedVertices();
        final_frag.insert(final_frag.end(), contract.begin(), contract.end());
      }
      for (CMGVertex v : overlap_vertices) {
        final_overlap.emplace_back(v.GetSource());
        auto con = v.GetContractedVertices();
        final_overlap.insert(final_overlap.end(), con.begin(), con.end());
      }
      Fragment f(MG, final_frag, final_overlap);
      if (std::find(_frags.at(m).begin(), _frags.at(m).end(), f)
          == _frags.at(m).end())
        _frags.at(m).emplace_back(f);
    }
    
    return _frags.at(m).size() - initial_count;
  }
  
}

