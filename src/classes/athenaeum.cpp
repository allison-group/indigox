#include <indigox/algorithm/graph/connectivity.hpp>
#include <indigox/algorithm/graph/paths.hpp>
#include <indigox/classes/angle.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/serialise.hpp>

#include <boost/dynamic_bitset.hpp>

#include <EASTL/iterator.h>
#include <EASTL/vector_set.h>
#include <algorithm>
#include <deque>
#include <fstream>
#include <iterator>
#include <numeric>
#include <vector>

namespace indigox {
  // ===========================================================================
  // == Fragment implementation ================================================
  // ===========================================================================

  struct Fragment::FragmentData {
    Molecule source_molecule;
    graph::CondensedMolecularGraph graph;
    std::vector<graph::CMGVertex> frag;
    std::vector<Fragment::OverlapVertex> overlap;
    std::vector<Fragment::AtmType> atoms;
    std::vector<Fragment::BndType> bonds;
    std::vector<Fragment::AngType> angles;
    std::vector<Fragment::DhdType> dihedrals;
    boost::dynamic_bitset<> supersets;
    boost::dynamic_bitset<> graph_mask;

    FragmentData() = default;

    template <class Archive> void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("mol", source_molecule),
              INDIGOX_SERIAL_NVP("graph", graph),
              INDIGOX_SERIAL_NVP("frag_verts", frag),
              INDIGOX_SERIAL_NVP("overlap_verts", overlap),
              INDIGOX_SERIAL_NVP("atoms", atoms),
              INDIGOX_SERIAL_NVP("bonds", bonds),
              INDIGOX_SERIAL_NVP("angles", angles),
              INDIGOX_SERIAL_NVP("dihedrals", dihedrals),
              INDIGOX_SERIAL_NVP("supersets", supersets),
              INDIGOX_SERIAL_NVP("mask", graph_mask));
    }
  };

  template <class Archive>
  void Fragment::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(Fragment);

  void _MGVertexToCMGVertex(std::vector<graph::MGVertex> &v_in,
                            std::vector<graph::CMGVertex> &v_out,
                            graph::CondensedMolecularGraph &g) {
    v_out.clear();
    eastl::vector_set<graph::MGVertex> contracted, relaxed;
    for (graph::MGVertex &v : v_in) {
      if (g.HasVertex(v))
        relaxed.emplace(v);
      else if (g.HasCondensedVertex(v))
        contracted.emplace(v);
      else
        throw std::runtime_error("Unknown vertex provided");
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

  Fragment::OverlapType _GetOverType(const graph::CMGVertex &v,
                                     graph::CondensedMolecularGraph &g,
                                     std::vector<graph::CMGVertex> &) {
    if (!g.HasVertex(v)) throw std::runtime_error("Issue with overlap type");
    return Fragment::OverlapType::GenericOverlap;
  }

  Fragment::Fragment(const Molecule &mol, std::vector<Atom> &frag,
                     std::vector<Atom> &overlap) {

    if (frag.empty()) { throw std::runtime_error("A fragment requires atoms"); }

    graph::MolecularGraph G = mol.GetGraph();
    std::vector<graph::MGVertex> frag_v, overlap_v;
    frag_v.reserve(frag.size());
    overlap_v.reserve(overlap.size());

    for (Atom atm : frag) frag_v.push_back(G.GetVertex(atm));
    for (Atom atm : overlap) overlap_v.push_back(G.GetVertex(atm));

    Fragment tmp(G, frag_v, overlap_v);
    m_data = tmp.m_data;
  }

  Fragment::Fragment(const graph::MolecularGraph &G,
                     std::vector<graph::MGVertex> &frag,
                     std::vector<graph::MGVertex> &overlap)
      : m_data(std::make_shared<FragmentData>()) {
    if (frag.empty()) throw std::runtime_error("A fragment needs vertices");

    m_data->source_molecule = G.GetMolecule();

    // Induce a new subgraph
    graph::MolecularGraph g = G;
    graph::CondensedMolecularGraph CG = g.GetCondensedGraph();
    std::vector<graph::CMGVertex> contracted_overlap, combined;
    _MGVertexToCMGVertex(frag, m_data->frag, CG);
    std::sort(m_data->frag.begin(), m_data->frag.end());
    _MGVertexToCMGVertex(overlap, contracted_overlap, CG);
    combined.assign(m_data->frag.begin(), m_data->frag.end());
    combined.insert(combined.end(), contracted_overlap.begin(),
                    contracted_overlap.end());
    m_data->graph = CG.Subgraph(combined);
    if (!m_data->graph.IsConnected())
      throw std::runtime_error("A fragment must be connected");
    m_data->atoms.assign(frag.begin(), frag.end());
    m_data->overlap.reserve(contracted_overlap.size());
    for (graph::CMGVertex &v : contracted_overlap)
      m_data->overlap.emplace_back(_GetOverType(v, m_data->graph, m_data->frag),
                                   v);
    std::sort(m_data->overlap.begin(), m_data->overlap.end());

    // Combine all MGVertices for checking purposes
    std::vector<AtmType> atm_check(frag.begin(), frag.end());
    atm_check.insert(atm_check.end(), overlap.begin(), overlap.end());
    long acceptable_pos = frag.size();

    // Get the molecule
    Molecule mol = frag.front().GetAtom().GetMolecule();

    // Set the graph_mask
    auto &all_v = CG.GetVertices();
    m_data->graph_mask = boost::dynamic_bitset<>(all_v.size());
    for (size_t i = 0; i < all_v.size(); ++i) {
      if (m_data->graph.HasVertex(all_v[i])) m_data->graph_mask.set(i);
    }

    // Get the bonds which are allowed
    std::deque<BndType> tmp_bnd;
    for (Bond bnd : mol.GetBonds()) {
      if (!bnd.HasType()) continue;
      AtmType a = G.GetVertex(bnd.GetAtoms()[0]);
      AtmType b = G.GetVertex(bnd.GetAtoms()[1]);
      auto a_pos = std::find(atm_check.begin(), atm_check.end(), a);
      if (a_pos == atm_check.end()) continue;
      auto b_pos = std::find(atm_check.begin(), atm_check.end(), b);
      if (b_pos == atm_check.end()) continue;
      size_t dangle = 0;
      if (std::distance(atm_check.begin(), a_pos) >= acceptable_pos) ++dangle;
      if (std::distance(atm_check.begin(), b_pos) >= acceptable_pos) ++dangle;
      // Bonds are allowed one dangling. Danglings go on end
      if (dangle > 1) continue;
      if (dangle)
        tmp_bnd.emplace_back(a, b);
      else
        tmp_bnd.emplace_front(a, b);
    }
    m_data->bonds.assign(tmp_bnd.begin(), tmp_bnd.end());

    // Get the angles which are allowed
    std::deque<AngType> tmp_ang;
    for (Angle ang : mol.GetAngles()) {
      if (!ang.HasType()) continue;
      AtmType a = G.GetVertex(ang.GetAtoms()[0]);
      AtmType b = G.GetVertex(ang.GetAtoms()[1]);
      AtmType c = G.GetVertex(ang.GetAtoms()[2]);
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
      if (dangle)
        tmp_ang.emplace_back(a, b, c);
      else
        tmp_ang.emplace_front(a, b, c);
    }
    m_data->angles.assign(tmp_ang.begin(), tmp_ang.end());

    // Get the dihedrals which are allowed
    std::deque<DhdType> tmp_dhd;
    for (Dihedral dhd : mol.GetDihedrals()) {
      if (!dhd.HasType()) continue;
      AtmType a = G.GetVertex(dhd.GetAtoms()[0]);
      AtmType b = G.GetVertex(dhd.GetAtoms()[1]);
      AtmType c = G.GetVertex(dhd.GetAtoms()[2]);
      AtmType d = G.GetVertex(dhd.GetAtoms()[3]);
      auto a_pos = std::find(atm_check.begin(), atm_check.end(), a);
      if (a_pos == atm_check.end()) continue;
      auto b_pos = std::find(atm_check.begin(), atm_check.end(), b);
      if (b_pos == atm_check.end()) continue;
      auto c_pos = std::find(atm_check.begin(), atm_check.end(), c);
      if (c_pos == atm_check.end()) continue;
      auto d_pos = std::find(atm_check.begin(), atm_check.end(), d);
      if (d_pos == atm_check.end()) continue;
      bool d1 = std::distance(atm_check.begin(), a_pos) < acceptable_pos;
      bool d2 = std::distance(atm_check.begin(), b_pos) < acceptable_pos;
      bool d3 = std::distance(atm_check.begin(), c_pos) < acceptable_pos;
      bool d4 = std::distance(atm_check.begin(), d_pos) < acceptable_pos;
      // Need two adjacent atoms to not be dangling
      if (!((d1 && d2) || (d2 && d3) || (d3 && d4))) continue;
      if (d1 || d2 || d3 || d4)
        tmp_dhd.emplace_front(a, b, c, d);
      else
        tmp_dhd.emplace_back(a, b, c, d);
    }
    m_data->dihedrals.assign(tmp_dhd.begin(), tmp_dhd.end());
  }

  const graph::CondensedMolecularGraph &Fragment::GetGraph() const {
    return m_data->graph;
  }

  const boost::dynamic_bitset<> &Fragment::GetSupersets() const {
    return m_data->supersets;
  }

  const std::vector<graph::CMGVertex> &Fragment::GetFragment() const {
    return m_data->frag;
  }

  size_t Fragment::Size() const { return m_data->frag.size(); }

  const std::vector<Fragment::OverlapVertex> &Fragment::GetOverlap() const {
    return m_data->overlap;
  }

  bool Fragment::IsFragmentVertex(const graph::CMGVertex &v) const {
    return (std::find(m_data->frag.begin(), m_data->frag.end(), v) !=
            m_data->frag.end());
  }

  bool Fragment::IsOverlapVertex(const graph::CMGVertex &v) const {
    auto pos = std::find_if(m_data->overlap.begin(), m_data->overlap.end(),
                            [&v](auto &u) { return u.second == v; });
    return pos != m_data->overlap.end();
  }

  const std::vector<Fragment::AtmType> &Fragment::GetAtoms() const {
    return m_data->atoms;
  }

  const std::vector<Fragment::BndType> &Fragment::GetBonds() const {
    return m_data->bonds;
  }

  const std::vector<Fragment::AngType> &Fragment::GetAngles() const {
    return m_data->angles;
  }

  const std::vector<Fragment::DhdType> &Fragment::GetDihedrals() const {
    return m_data->dihedrals;
  }

  bool Fragment::operator==(const Fragment &frag) const {
    if (m_data == frag.m_data) return true;
    if (m_data->frag != frag.m_data->frag) return false;
    if (m_data->overlap != frag.m_data->overlap) return false;
    return true;
  }

  bool Fragment::operator<(const Fragment &frag) const {
    if (m_data == frag.m_data) return false;
    if (m_data->frag < frag.m_data->frag) return true;
    return m_data->overlap < frag.m_data->overlap;
  }

  bool Fragment::operator>(const Fragment &frag) const {
    if (m_data == frag.m_data) return false;
    if (m_data->frag > frag.m_data->frag) return true;
    return m_data->overlap > frag.m_data->overlap;
  }

  // ===========================================================================
  // == Athenaeum implementation ===============================================
  // ===========================================================================

  using AthSettings = Athenaeum::Settings;

  struct Athenaeum::Impl {

    std::bitset<(uint8_t)Settings::BoolCount> bool_parameters;
    std::array<int32_t,
               (uint8_t)Settings::IntCount - (uint8_t)Settings::BoolCount - 1>
        int_parameters;

    Forcefield ff;
    MoleculeFragments fragments;

    Impl() = default;
    Impl(const Forcefield &f) : bool_parameters(0), ff(f) {}

    template <class Archive> void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("bool_settings", bool_parameters),
              INDIGOX_SERIAL_NVP("int_settings", int_parameters),
              INDIGOX_SERIAL_NVP("forcefield", ff),
              INDIGOX_SERIAL_NVP("fragments", fragments));
    }
  };

  // Default settings

  void Athenaeum::DefaultSettings() {
    SetInt(AthSettings::MoleculeSizeLimit, 40);
    SetInt(AthSettings::MinimumFragmentSize, 1);
    SetInt(AthSettings::MaximumFragmentSize, 40);
    SetInt(AthSettings::OverlapLength, 2);
    SetInt(AthSettings::CycleSize, 8);
  }

  bool Athenaeum::GetBool(AthSettings param) {
    if (param >= AthSettings::BoolCount)
      throw std::runtime_error("Not a boolean parameter");
    return m_data->bool_parameters.test((uint8_t)param);
  }

  void Athenaeum::SetBool(AthSettings param) {
    if (param >= AthSettings::BoolCount)
      throw std::runtime_error("Not a boolean parameter");
    m_data->bool_parameters.set((uint8_t)param);
  }

  void Athenaeum::UnsetBool(AthSettings param) {
    if (param >= AthSettings::BoolCount)
      throw std::runtime_error("Not a boolean parameter");
    m_data->bool_parameters.reset((uint8_t)param);
  }

  int32_t Athenaeum::GetInt(AthSettings param) {
    uint8_t offset = 1 + (uint8_t)AthSettings::BoolCount;
    if (param <= AthSettings::BoolCount || param >= AthSettings::IntCount)
      throw std::runtime_error("Not an integer parameter");
    return m_data->int_parameters[(uint8_t)param - offset];
  }

  void Athenaeum::SetInt(AthSettings param, int32_t value) {
    uint8_t offset = 1 + (uint8_t)AthSettings::BoolCount;
    if (param <= AthSettings::BoolCount || param >= AthSettings::IntCount)
      throw std::runtime_error("Not an integer parameter");
    m_data->int_parameters[(uint8_t)param - offset] = value;
  }

  Athenaeum::Athenaeum(const Forcefield &ff)
      : m_data(std::make_shared<Impl>(ff)) {
    DefaultSettings();
  }

  Athenaeum::Athenaeum(const Forcefield &ff, int32_t overlap) : Athenaeum(ff) {
    SetInt(AthSettings::OverlapLength, overlap);
  }

  // cycle overlap currently ignored
  Athenaeum::Athenaeum(const Forcefield &ff, int32_t overlap, int32_t)
      : Athenaeum(ff) {
    SetInt(AthSettings::OverlapLength, overlap);
  }

  template <class Archive>
  void Athenaeum::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(Athenaeum);

  size_t Athenaeum::NumFragments() const {
    return std::accumulate(
        m_data->fragments.begin(), m_data->fragments.end(), 0,
        [](size_t i, auto &j) { return i + j.second.size(); });
  }

  size_t Athenaeum::NumFragments(const Molecule &mol) const {
    auto pos = m_data->fragments.find(mol);
    if (pos == m_data->fragments.end()) return 0;
    return pos->second.size();
  }

  const Athenaeum::MoleculeFragments &Athenaeum::GetFragments() const {
    return m_data->fragments;
  }

  const Athenaeum::FragContain &
  Athenaeum::GetFragments(const Molecule &mol) const {
    auto pos = m_data->fragments.find(mol);
    if (pos == m_data->fragments.end())
      throw std::runtime_error("No fragmenst for molecule available");
    return pos->second;
  }

  bool Athenaeum::HasFragments(const Molecule &mol) const {
    return m_data->fragments.find(mol) != m_data->fragments.end();
  }

  const Forcefield &Athenaeum::GetForcefield() const { return m_data->ff; }

  void Athenaeum::SortAndMask(const Molecule &mol) {
    auto pos = m_data->fragments.find(mol);
    FragContain &frags = pos->second;

    // Sort based on size
    std::sort(frags.begin(), frags.end(), [](Fragment &a, Fragment &b) {
      return a.GetGraph().NumVertices() < b.GetGraph().NumVertices();
    });

    for (size_t i = 0; i < frags.size(); ++i) {
      Fragment f = frags[i];
      f.m_data->supersets = boost::dynamic_bitset<>(frags.size());
      for (size_t j = i + 1; j < frags.size(); ++j) {
        if (f.m_data->graph_mask.is_proper_subset_of(
                frags[j].m_data->graph_mask))
          f.m_data->supersets.set(j);
      }
    }
  }

  bool Athenaeum::AddFragment(const Fragment &frag) {
    // Check that the fragment matches the molecule
    Molecule mol = frag.m_data->source_molecule;
    graph::MolecularGraph MG = mol.GetGraph();
    graph::CondensedMolecularGraph CG = MG.GetCondensedGraph();
    graph::CondensedMolecularGraph fg = frag.GetGraph();
    while (fg.IsSubgraph()) { fg = fg.GetSuperGraph(); }
    if (CG != fg) return false;

    // Check that the molecule forcefield matchs the athenaeum forcefield
    if (!mol.HasForcefield()) return false;
    if (mol.GetForcefield() != m_data->ff) return false;

    auto pos = m_data->fragments.emplace(mol, FragContain());
    pos.first->second.emplace_back(frag);
    SortAndMask(mol);
    return true;
  }

  bool CanCutEdge(graph::CMGEdge &e, graph::CondensedMolecularGraph &g) {
    if (!g.HasEdge(e)) return false;
    // Only single or aromatic
    Bond bnd = e.GetSource().GetBond();
    if (bnd.GetOrder() != BondOrder::SINGLE &&
        bnd.GetOrder() != BondOrder::AROMATIC)
      return false;
    // Is a hetero bond
    if (bnd.GetAtoms()[0].GetElement() != "C" &&
        bnd.GetAtoms()[1].GetElement() != "C")
      return false;
    // cyclisation rules
    return true;
  }

  size_t Athenaeum::AddAllFragments(const Molecule &mol) {
    using namespace indigox::graph;
    // Perform checks
    if (!mol.HasForcefield())
      throw std::runtime_error(
          "Attempting to fragment unparameterised molecule");
    if (mol.GetForcefield() != GetForcefield()) {
      std::stringstream message;
      message << "Forcefield mismatch. Expected " << GetForcefield().GetName()
              << " but molecule uses " << mol.GetForcefield().GetName();
      throw std::runtime_error(message.str());
    }

    int32_t mol_size_lim = GetInt(AthSettings::MoleculeSizeLimit);
    if (mol.NumAtoms() > mol_size_lim) {
      std::stringstream message;
      message << "Molecule too large to automagically fragment. Mol atoms: " << mol.NumAtoms()
              << ", Max atoms allowed: " << mol_size_lim;
      throw std::runtime_error(message.str());
    }

    auto pos = m_data->fragments.emplace(mol, FragContain());
    size_t initial_count = pos.first->second.size();

    // Get all the subgraphs of the molecule's condensed graph
    MolecularGraph MG = mol.GetGraph();
    CondensedMolecularGraph CG = MG.GetCondensedGraph();
    eastl::vector_set<CMGVertex> all_vertices(CG.GetVertices().begin(),
                                              CG.GetVertices().end());
    eastl::vector_set<CMGEdge> all_edges(CG.GetEdges().begin(),
                                         CG.GetEdges().end());

    // Decide if each subgraph can be made into a fragment
    CondensedMolecularGraph sub;
    algorithm::ConnectedSubgraphs subgraph_generator(CG);
    while (subgraph_generator(sub)) {
      // Sort the vertices/edges of CG into not in sub and in sub
      eastl::vector_set<CMGVertex> sub_vertices(sub.GetVertices().begin(),
                                                sub.GetVertices().end());
      eastl::vector_set<CMGEdge> sub_edges(sub.GetEdges().begin(),
                                           sub.GetEdges().end());
      std::vector<CMGVertex> other_vertices;
      std::vector<CMGEdge> other_edges;
      std::set_difference(all_vertices.begin(), all_vertices.end(),
                          sub_vertices.begin(), sub_vertices.end(),
                          std::back_inserter(other_vertices));
      std::set_difference(all_edges.begin(), all_edges.end(), sub_edges.begin(),
                          sub_edges.end(), std::back_inserter(other_edges));

      // Determine which edges are cut
      // cut_edge.second is true if source vertex in other_vertices, false if
      // not
      std::vector<std::pair<CMGEdge, bool>> cut_edges;
      for (CMGEdge e : other_edges) {
        CMGVertex u = CG.GetSourceVertex(e);
        CMGVertex v = CG.GetTargetVertex(e);
        bool has_u = sub.HasVertex(u);
        bool has_v = sub.HasVertex(v);
        if (has_u && has_v)
          throw std::runtime_error("WTF?!");
        else if (has_u && !has_v)
          cut_edges.emplace_back(e, false);
        else if (!has_u && has_v)
          cut_edges.emplace_back(e, true);
        else
          continue;
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
          if (!path.empty() &&
              (int32_t)path.size() <= GetInt(AthSettings::OverlapLength))
            overlap_vertices.emplace(u);
        }
      }

      // every leaf in overlap must have minimum path length of _overlap to
      // each vertex in fragment
      std::vector<CMGVertex> fragoververt(sub_vertices.begin(),
                                          sub_vertices.end());
      fragoververt.insert(fragoververt.end(), overlap_vertices.begin(),
                          overlap_vertices.end());
      CondensedMolecularGraph withoverlap = CG.Subgraph(fragoververt);
      CondensedMolecularGraph::ComponentContain tmp;
      if (algorithm::ConnectedComponents(withoverlap, tmp) > 1) continue;
      bool bad_overlaps = false;
      for (CMGVertex u : overlap_vertices) {
        if (withoverlap.Degree(u) > 1) continue;
        for (CMGVertex v : sub_vertices) {
          auto path = algorithm::ShortestPath(withoverlap, u, v);
          if ((int32_t)path.size() < GetInt(AthSettings::OverlapLength))
            bad_overlaps = true;
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
      if (std::find(pos.first->second.begin(), pos.first->second.end(), f) ==
          pos.first->second.end())
        pos.first->second.emplace_back(f);
    }
    if (pos.first->second.size() != initial_count) { SortAndMask(mol); }
    return pos.first->second.size() - initial_count;
  }

  void SaveAthenaeum(const Athenaeum &a, const std::string& path) {
    using Archive = cereal::PortableBinaryOutputArchive;
    std::ofstream os(path);
    if (!os.is_open()) throw std::runtime_error("Unable to open output stream");
    Archive archive(os);
    std::string stype("Athenaeum");
    archive(stype, a);
  }

  Athenaeum LoadAthenaeum(std::string path) {
    using Archive = cereal::PortableBinaryInputArchive;
    std::ifstream is(path);
    if (!is.is_open()) throw std::runtime_error("Unable to open input stream");
    std::string stype;
    Archive archive(is);
    archive(stype);
    if (stype != "Athenaeum") throw std::runtime_error("Not an Athenaeum file");
    Athenaeum a;
    archive(a);
    return a;
  }

  bool Athenaeum::operator==(const Athenaeum &a) const {
    return m_data == a.m_data;
  }

  bool Athenaeum::operator<(const Athenaeum &ath) const {
    return m_data < ath.m_data;
  }

  bool Athenaeum::operator>(const Athenaeum &ath) const {
    return m_data > ath.m_data;
  }

  std::ostream &operator<<(std::ostream &os, const Fragment &frag) {
    if (frag) {
      os << "Fragment(" << frag.m_data->frag.size() << " core, "
         << frag.m_data->overlap.size() << " overlap)";
    }
    return os;
  }

} // namespace indigox
