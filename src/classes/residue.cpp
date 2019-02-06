#include <indigox/classes/atom.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/molecule_impl.hpp>
#include <indigox/classes/residue.hpp>
#include <indigox/graph/molecular.hpp>

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

  void Residue::Impl::DetermineType() { type = ResidueType::NonSpecific; }

  ResidueType Residue::GetType() {
    _sanity_check_(*this);
    if (m_data->type == ResidueType::Unspecified) { m_data->DetermineType(); }
    return m_data->type;
  }

  bool Residue::IsAminoAcid() {
    _sanity_check_(*this);
    return GetType() == ResidueType::AminoAcid;
  }

  bool Residue::IsAlphaAminoAcid() {
    _sanity_check_(*this);
    return false;
  }

  bool Residue::IsBetaAminoAcid() {
    _sanity_check_(*this);
    return false;
  }

  bool Residue::IsGammaAminoAcid() {
    _sanity_check_(*this);
    return false;
  }

  bool Residue::IsDeltaAminoAcid() {
    _sanity_check_(*this);
    return false;
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

} // namespace indigox
