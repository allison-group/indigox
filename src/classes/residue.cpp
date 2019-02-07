#include <indigox/classes/atom.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
#include <indigox/classes/residue.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/algorithm/graph/paths.hpp>

#ifndef INDIGOX_DISABLE_SANITY_CHECKS
#define _sanity_check_(x)                                                      \
  if (!x)                                                                      \
  throw std::runtime_error(                                                    \
      "Attempting to access data from invalid residue instance")
#else
#define _sanity_check_(x)
#endif

namespace indigox {

  // =======================================================================
  // == RESIDUE ============================================================
  // =======================================================================

  Residue::Impl::Impl(const std::vector<Atom> &atms, const Molecule &mol)
      : type(ResidueType::Unspecified), molecule(mol) {
    atoms.insert(atms.begin(), atms.end());
    graph::MolecularGraph g = mol.GetGraph();
    std::vector<graph::MGVertex> vertices;
    vertices.reserve(atms.size());
    for (Atom atm : atms) vertices.push_back(g.GetVertex(atm));
    residue_graph = g.Subgraph(vertices);
  }

  Residue::Residue(const std::vector<Atom> &atms, const Molecule &mol)
      : m_data(std::make_shared<Impl>(atms, mol)) {}

  template <typename Archive>
  void Residue::Impl::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("type", type),
            INDIGOX_SERIAL_NVP("atoms", atoms),
            INDIGOX_SERIAL_NVP("molecule", molecule),
            INDIGOX_SERIAL_NVP("graph", residue_graph));
  }

  template <typename Archive>
  void Residue::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }

  bool Residue::HasAtom(const Atom &atom) const {
    _sanity_check_(*this);
    return m_data->atoms.find(atom) != m_data->atoms.end();
  }

  void Residue::Impl::DetermineType() {
    if (AminoAcidTest()) type = ResidueType::AminoAcid;
    else type = ResidueType::NonSpecific;
  }

  ResidueType Residue::GetType() {
    _sanity_check_(*this);
    if (m_data->type == ResidueType::Unspecified) { m_data->DetermineType(); }
    return m_data->type;
  }

  bool Residue::IsAminoAcid() {
    _sanity_check_(*this);
    return m_data->AminoAcidTest();
  }

  bool Residue::IsAlphaAminoAcid() {
    _sanity_check_(*this);
    if (!IsAminoAcid()) return false;
    return m_data->cache_aa_length.find(1) != m_data->cache_aa_length.end();
  }

  bool Residue::IsBetaAminoAcid() {
    _sanity_check_(*this);
    if (!IsAminoAcid()) return false;
    return m_data->cache_aa_length.find(2) != m_data->cache_aa_length.end();
  }

  bool Residue::IsGammaAminoAcid() {
    _sanity_check_(*this);
    if (!IsAminoAcid()) return false;
    return m_data->cache_aa_length.find(3) != m_data->cache_aa_length.end();
  }

  bool Residue::IsDeltaAminoAcid() {
    _sanity_check_(*this);
    if (!IsAminoAcid()) return false;
    return m_data->cache_aa_length.find(4) != m_data->cache_aa_length.end();
  }

  bool Residue::IsSugar() {
    _sanity_check_(*this);
    return GetType() == ResidueType::Sugar;
  }

  bool Residue::IsLipid() {
    _sanity_check_(*this);
    return GetType() == ResidueType::Lipid;
  }

  bool Residue::IsNucleicAcid() {
    _sanity_check_(*this);
    return GetType() == ResidueType::NucleicAcid;
  }

  const Residue::ResidueAtoms &Residue::GetAtoms() const {
    _sanity_check_(*this);
    return m_data->atoms;
  }

  bool Residue::Impl::AminoAcidTest() {
    if (!cache_aa_length.empty()) {
      return cache_aa_length.size() > 1;
    }
    cache_aa_length.insert(-1);
    
    graph::MolecularGraph G = molecule.GetGraph();
    std::vector<graph::MGVertex> carbons, nitrogens;
    for (Atom atm : atoms) {
      graph::MGVertex v = residue_graph.GetVertex(atm);
      if (atm.GetElement() == "C") {
        bool carbonyl = false, oxygen = false;
        for (graph::MGVertex u : residue_graph.GetNeighbours(v)) {
          Bond bnd = residue_graph.GetEdge(u, v).GetBond();
          if (!carbonyl && bnd.IsCarbonylBond()) carbonyl = true;
          else if (!oxygen && u.GetAtom().GetElement() == "O") oxygen = true;
        }
        if (carbonyl && (oxygen || residue_graph.Degree(v) + 1 == G.Degree(v))) {
          carbons.push_back(v);
        }
      }
      if (atm.GetElement() == "N") {
        if (G.Degree(v) - 1 == residue_graph.Degree(v)) nitrogens.push_back(v);
        else {
          for (graph::MGVertex u : residue_graph.GetNeighbours(v)) {
            if (u.GetAtom().GetElement() == "H") {
              nitrogens.push_back(v);
              break;
            }
          }
        }
      }
    }

    if (carbons.empty() || nitrogens.empty()) return false;
    
    for (graph::MGVertex source : carbons) {
      for (graph::MGVertex target : nitrogens) {
        auto path = algorithm::ShortestPath(residue_graph, source, target);
        std::vector<graph::MGVertex> v_path;
        for (graph::MGEdge edge : path) {
          graph::MGVertex begin = residue_graph.GetSourceVertex(edge);
          graph::MGVertex end = residue_graph.GetTargetVertex(edge);
          if (v_path.empty() && begin == source) v_path.emplace_back(end);
          else if (v_path.empty() && end == source) v_path.emplace_back(begin);
          else if (begin == v_path.back() && end != target) v_path.emplace_back(end);
          else if (end == v_path.back() && begin != target) v_path.emplace_back(begin);
        }
        bool is_amino_acid_path = !v_path.empty();
        for (graph::MGVertex v : v_path){
          if (v.GetAtom().GetElement() != "C") is_amino_acid_path = false;
        }
        if (is_amino_acid_path) cache_aa_length.insert(v_path.size());
      }
    }
    
    return cache_aa_length.size() > 1;
  }

} // namespace indigox
