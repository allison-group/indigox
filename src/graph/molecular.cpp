#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/serialise.hpp>

#include <algorithm>
#include <limits>

namespace indigox::graph {

  struct MGVertex::MGVertexData {
    Atom atom;
    MolecularGraph graph;

    MGVertexData() = default;
    MGVertexData(const Atom &a, const MolecularGraph &g) : atom(a), graph(g) {}

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("atom", atom),
              INDIGOX_SERIAL_NVP("graph", graph));
    }
  };

  struct MGEdge::MGEdgeData {
    Bond bond;
    MolecularGraph graph;

    MGEdgeData() = default;
    MGEdgeData(const Bond &b, const MolecularGraph &g) : bond(b), graph(g) {}

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("bond", bond),
              INDIGOX_SERIAL_NVP("graph", graph));
    }
  };

  struct MolecularGraph::Impl {
    AtomMap atom_vertices;
    BondMap bond_edges;
    Molecule molecule;
    CondensedMolecularGraph condensed_graph;
    MolecularGraph super_graph;

    Impl() = default;
    Impl(const Molecule &mol) : molecule(mol) {}

    template <typename Archive>
    void serialise(Archive &archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("atom_map", atom_vertices),
              INDIGOX_SERIAL_NVP("bond_map", bond_edges),
              INDIGOX_SERIAL_NVP("molecule", molecule),
              INDIGOX_SERIAL_NVP("condensed_graph", condensed_graph),
              INDIGOX_SERIAL_NVP("supergraph", super_graph));
    }
  };

  // ==================================================================
  // == SERIALISATION =================================================
  // ==================================================================

  template <typename Archive>
  void MGVertex::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(MGVertex);

  template <typename Archive>
  void MGEdge::serialise(Archive &archive, const uint32_t) {
    archive(INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(MGEdge);

  template <typename Archive>
  void MolecularGraph::serialise(Archive &archive, const uint32_t) {
    archive(
        INDIGOX_SERIAL_NVP("base_graph", cereal::base_class<graph_type>(this)),
        INDIGOX_SERIAL_NVP("data", m_data));
  }
  INDIGOX_SERIALISE(MolecularGraph);

  // ====================================================================
  // == CONSTRUCTION ====================================================
  // ====================================================================

  MGVertex::MGVertex(const Atom &a, const MolecularGraph &graph)
      : m_data(std::make_shared<MGVertexData>(a, graph)) {}

  MGEdge::MGEdge(const Bond &b, const MolecularGraph &graph)
      : m_data(std::make_shared<MGEdgeData>(b, graph)) {}

  MolecularGraph::MolecularGraph(const Molecule &mol)
      : BaseGraph<MGVertex, MGEdge, MolecularGraph>(),
        m_data(std::make_shared<Impl>(mol)) {}

  // =====================================================================
  // == GETTING ==========================================================
  // =====================================================================

  const Atom &MGVertex::GetAtom() const { return m_data->atom; }
  const MolecularGraph &MGVertex::GetGraph() const { return m_data->graph; }

  const Bond &MGEdge::GetBond() const { return m_data->bond; }
  const MolecularGraph &MGEdge::GetGraph() const { return m_data->graph; }

  const MGEdge &MolecularGraph::GetEdge(const Bond &bnd) const {
    if (!HasEdge(bnd)) throw std::out_of_range("No such edge");
    return m_data->bond_edges.at(bnd);
  }

  const MGVertex &MolecularGraph::GetVertex(const Atom &atm) const {
    if (!HasVertex(atm)) throw std::out_of_range("No such vertex");
    return m_data->atom_vertices.at(atm);
  }

  const Molecule &MolecularGraph::GetMolecule() const {
    if (IsSubgraph())
      throw std::runtime_error("Subgraphs do not relate to a molecule");
    return m_data->molecule;
  }

  const MolecularGraph &MolecularGraph::GetSuperGraph() const {
    if (!IsSubgraph())
      throw std::runtime_error("Cannot get supergraph as not a subgraph");
    return m_data->super_graph;
  }

  bool MolecularGraph::IsSubgraph() const { return bool(m_data->super_graph); }

  // ===================================================================
  // == Subgraph Generation ============================================
  // ===================================================================
  MolecularGraph MolecularGraph::Subgraph(std::vector<MGVertex> &verts) {
    MolecularGraph G;
    G.m_data = std::make_shared<Impl>();
    G.m_data->super_graph = *this;

    for (const MGVertex &v : verts) {
      if (!HasVertex(v)) throw std::runtime_error("Non-member vertex");
      G.m_data->atom_vertices.emplace(v.GetAtom(), v);
      G.graph_type::AddVertex(v);
    }

    for (const MGEdge &e : GetEdges()) {
      MGVertex u = GetSourceVertex(e);
      MGVertex v = GetTargetVertex(e);
      if (!G.HasVertex(u) || !G.HasVertex(v)) { continue; }
      G.m_data->bond_edges.emplace(e.GetBond(), e);
      G.graph_type::AddEdge(u, v, e);
    }
    return G;
  }

  MolecularGraph MolecularGraph::Subgraph(std::vector<MGVertex> &verts,
                                          std::vector<MGEdge> &edges) {
    MolecularGraph G;
    G.m_data = std::make_shared<Impl>();
    G.m_data->super_graph = *this;

    for (const MGVertex &v : verts) {
      if (!HasVertex(v)) throw std::runtime_error("Non-member vertex");
      G.m_data->atom_vertices.emplace(v.GetAtom(), v);
      G.graph_type::AddVertex(v);
    }

    for (const MGEdge &e : edges) {
      if (!HasEdge(e)) throw std::runtime_error("Non-member edge");
      MGVertex u = GetSourceVertex(e);
      MGVertex v = GetTargetVertex(e);
      if (!G.HasVertex(u) || !G.HasVertex(v)) continue;
      G.m_data->bond_edges.emplace(e.GetBond(), e);
      G.graph_type::AddEdge(u, v, e);
    }
    return G;
  }

  // =====================================================================
  // == STRUCTURE MODIFICATION ===========================================
  // =====================================================================

  MGVertex MolecularGraph::AddVertex(const Atom &atm) {
    MGVertex v(atm, *this);
    m_data->atom_vertices.emplace(atm, v);
    graph_type::AddVertex(v);
    return v;
  }

  void MolecularGraph::RemoveVertex(const MGVertex &v) {
    m_data->atom_vertices.erase(v.GetAtom());
    for (const MGVertex &u : graph_type::GetNeighbours(v)) {
      MGEdge e = GetEdge(u, v);
      m_data->bond_edges.erase(e.GetBond());
    }
    graph_type::RemoveVertex(v);
  }

  MGEdge MolecularGraph::AddEdge(const Bond &bnd) {
    Bond::BondAtoms atoms = bnd.GetAtoms();
    // Edges will only be added when a bond is added to the molecule
    MGVertex u = GetVertex(atoms[0]);
    MGVertex v = GetVertex(atoms[1]);
    MGEdge e(bnd, *this);
    m_data->bond_edges.emplace(bnd, e);
    graph_type::AddEdge(u, v, e);
    return e;
  }

  void MolecularGraph::RemoveEdge(const MGEdge &e) {
    m_data->bond_edges.erase(e.GetBond());
    graph_type::RemoveEdge(e);
  }

  void MolecularGraph::RemoveEdge(const MGVertex &u, const MGVertex &v) {
    MGEdge e = GetEdge(u, v);
    RemoveEdge(e);
  }

  bool MolecularGraph::HasVertex(const Atom &v) const {
    return m_data->atom_vertices.find(v) != m_data->atom_vertices.end();
  }

  bool MolecularGraph::HasEdge(const Bond &e) const {
    return m_data->bond_edges.find(e) != m_data->bond_edges.end();
  }

  const CondensedMolecularGraph &MolecularGraph::GetCondensedGraph() {
    if (!m_data->molecule.IsFrozen())
      throw std::runtime_error("Can only condense a frozen molecule");
    if (!m_data->condensed_graph) m_data->condensed_graph = Condense(*this);
    return m_data->condensed_graph;
  }

  // =====================================================================
  // == OPERATORS ========================================================
  // =====================================================================

  bool MGVertex::operator==(const MGVertex &v) const {
    return m_data->atom == v.m_data->atom;
  }

  bool MGVertex::operator<(const MGVertex &v) const {
    return m_data->atom < v.m_data->atom;
  }

  bool MGVertex::operator>(const MGVertex &v) const {
    return m_data->atom > v.m_data->atom;
  }

  std::ostream &operator<<(std::ostream &os, const MGVertex &v) {
    if (v) { os << "MGVertex(" << v.GetAtom().GetIndex() + 1 << ")"; }
    return os;
  }

  bool MGEdge::operator==(const MGEdge &v) const {
    return m_data->bond == v.m_data->bond;
  }

  bool MGEdge::operator<(const MGEdge &v) const {
    return m_data->bond < v.m_data->bond;
  }

  bool MGEdge::operator>(const MGEdge &v) const {
    return m_data->bond > v.m_data->bond;
  }

  std::ostream &operator<<(std::ostream &os, const MGEdge &e) {
    if (e) { os << "MGEdge(" << e.GetBond().GetIndex() + 1 << ")"; }
    return os;
  }

  bool MolecularGraph::operator==(const MolecularGraph &g) const {
    return m_data == g.m_data;
  }

  bool MolecularGraph::operator<(const MolecularGraph &g) const {
    return m_data < g.m_data;
  }

  bool MolecularGraph::operator>(const MolecularGraph &g) const {
    return m_data > g.m_data;
  }

  std::ostream &operator<<(std::ostream &os, const MolecularGraph &G) {
    if (G) {
      os << "MolecularGraph(" << G.NumVertices() << " vertices, "
         << G.NumEdges() << " edges)";
    }
    return os;
  }

} // namespace indigox::graph
