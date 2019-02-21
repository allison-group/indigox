#include <indigox/algorithm/cherrypicker.hpp>
#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/classes/angle.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/combinatronics.hpp>

#include <rilib/RI.h>

#include <boost/dynamic_bitset.hpp>

#include <algorithm>
#include <map>
#include <vector>

namespace indigox::algorithm {

  struct Uint64AttrComparator : public rilib::AttributeComparator {
    Uint64AttrComparator(){};
    virtual bool compare(void *attr1, void *attr2) {
      return (*((uint64_t *)attr1) - *((uint64_t *)attr2)) == 0;
    }
    virtual int compareint(void *attr1, void *attr2) {
      return *((uint64_t *)attr1) - *((uint64_t *)attr2);
    }
  };

  struct Uint32AttrComparator : public rilib::AttributeComparator {
    Uint32AttrComparator(){};
    virtual bool compare(void *attr1, void *attr2) {
      return (*((uint32_t *)attr1) - *((uint32_t *)attr2)) == 0;
    }
    virtual int compareint(void *attr1, void *attr2) {
      return *((uint32_t *)attr1) - *((uint32_t *)attr2);
    }
  };

  std::unique_ptr<rilib::Graph>
  CMGToRIGraph(graph::CondensedMolecularGraph &cmg, graph::EdgeIsoMask emask,
               graph::VertexIsoMask vmask) {
    std::unique_ptr<rilib::Graph> G = std::make_unique<rilib::Graph>();

    std::vector<graph::CMGVertex> cmg_v = cmg.GetVertices();

    // Build the vertices
    G->nof_nodes = cmg_v.size();
    G->nodes_attrs = (void **)malloc(G->nof_nodes * sizeof(void *));
    G->out_adj_sizes = (int *)calloc(G->nof_nodes, sizeof(int));
    G->in_adj_sizes = (int *)calloc(G->nof_nodes, sizeof(int));
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->nodes_attrs[i] = (uint64_t *)malloc(sizeof(uint64_t));
      *((uint64_t *)G->nodes_attrs[i]) =
          (cmg_v[i].GetIsomorphismMask() & vmask).to_uint64();
      G->out_adj_sizes[i] += cmg.Degree(cmg_v[i]);
      G->in_adj_sizes[i] += cmg.Degree(cmg_v[i]);
    }

    // Build the edges
    G->out_adj_list = (int **)malloc(G->nof_nodes * sizeof(int *));
    G->in_adj_list = (int **)malloc(G->nof_nodes * sizeof(int *));
    G->out_adj_attrs = (void ***)malloc(G->nof_nodes * sizeof(void **));

    int *ink = (int *)calloc(G->nof_nodes, sizeof(int));
    for (int i = 0; i < G->nof_nodes; ++i) {
      G->in_adj_list[i] = (int *)calloc(G->in_adj_sizes[i], sizeof(int));
    }

    for (int i = 0; i < G->nof_nodes; ++i) {
      G->out_adj_list[i] = (int *)calloc(G->in_adj_sizes[i], sizeof(int));
      G->out_adj_attrs[i] =
          (void **)malloc(G->out_adj_sizes[i] * sizeof(void *));
      std::vector<graph::CMGVertex> nbrs = cmg.GetNeighbours(cmg_v[i]);
      for (int j = 0; j < G->out_adj_sizes[i]; ++j) {
        int idx = std::distance(cmg_v.begin(),
                                std::find(cmg_v.begin(), cmg_v.end(), nbrs[j]));
        G->out_adj_list[i][j] = idx;
        G->out_adj_attrs[i][j] = (uint32_t *)malloc(sizeof(uint32_t));
        *((uint32_t *)G->out_adj_attrs[i][j]) =
            (cmg.GetEdge(cmg_v[i], nbrs[j]).GetIsomorphismMask() & emask)
                .to_uint32();
        G->in_adj_list[idx][ink[idx]] = i;
        ink[idx]++;
      }
    }

    free(ink);
    return G;
  }

  using CPSet = CherryPicker::Settings;

  CherryPicker::CherryPicker(Forcefield &ff) : _ff(ff), bool_parameters(0) {
    DefaultSettings();
  }

  void CherryPicker::DefaultSettings() {
    SetBool(CPSet::VertexElement);
    SetBool(CPSet::VertexFormalCharge);
    SetBool(CPSet::VertexCondensed);
    SetBool(CPSet::VertexDegree);

    SetBool(CPSet::EdgeBondOrder);
    SetBool(CPSet::EdgeDegree);

    SetBool(CPSet::AllowDanglingBonds);
    SetBool(CPSet::AllowDanglingAngles);
    SetBool(CPSet::AllowDanglingDihedrals);
    SetBool(CPSet::UseRISubgraphMatching);

    SetInt(CPSet::MinimumFragmentSize, 4);
    SetInt(CPSet::MaximumFragmentSize, -1);
  }

  bool CherryPicker::GetBool(CPSet param) {
    if (param >= CPSet::BoolCount)
      throw std::runtime_error("Not a boolean parameter");
    return bool_parameters.test((uint8_t)param);
  }

  void CherryPicker::SetBool(CPSet param) {
    if (param >= CPSet::BoolCount)
      throw std::runtime_error("Not a boolean parameter");
    bool_parameters.set((uint8_t)param);
  }

  void CherryPicker::UnsetBool(CPSet param) {
    if (param >= CPSet::BoolCount)
      throw std::runtime_error("Not a boolean parameter");
    bool_parameters.reset((uint8_t)param);
  }

  int32_t CherryPicker::GetInt(CPSet param) {
    uint8_t offset = 1 + (uint8_t)CPSet::BoolCount;
    if (param <= CPSet::BoolCount || param >= CPSet::IntCount)
      throw std::runtime_error("Not an integer parameter");
    return int_parameters[(uint8_t)param - offset];
  }

  void CherryPicker::SetInt(CPSet param, int32_t value) {
    uint8_t offset = 1 + (uint8_t)CPSet::BoolCount;
    if (param <= CPSet::BoolCount || param >= CPSet::IntCount)
      throw std::runtime_error("Not an integer parameter");
    int_parameters[(uint8_t)param - offset] = value;
  }

  bool CherryPicker::AddAthenaeum(Athenaeum &library) {
    if (library.GetForcefield() != _ff) return false;
    _libs.push_back(library);
    return true;
  }

  bool CherryPicker::RemoveAthenaeum(Athenaeum &library) {
    auto pos = std::find(_libs.begin(), _libs.end(), library);
    if (pos != _libs.end()) _libs.erase(pos);
    return pos != _libs.end();
  }

  using CMGV = graph::CMGVertex;
  using CMGE = graph::CMGEdge;
  using CMGS = graph::CondensedMolecularGraph;
  using U = graph::Undirected;
  using GL = graph::GraphLabel;

  struct CherryPickerCallback : public CMGCallback {
    using GraphType = graph::CondensedMolecularGraph;
    using BaseType = CMGCallback;
    using CorrespondenceMap = BaseType::CorrespondenceMap;

    using VertMasks = eastl::vector_map<CMGV, graph::VertexIsoMask>;
    using EdgeMasks = eastl::vector_map<CMGE, graph::EdgeIsoMask>;

    CherryPicker &cherrypicker;
    GraphType small;
    GraphType large;
    VertMasks &vmasks_large;
    EdgeMasks &emasks_large;
    VertMasks vmasks_small;
    EdgeMasks emasks_small;
    ParamMolecule pmol;
    Fragment frag;
    bool has_mapping;

    CherryPickerCallback(CherryPicker &cp, GraphType &l, VertMasks &vl,
                         EdgeMasks &el, ParamMolecule &p, Fragment &f,
                         graph::VertexIsoMask vertmask,
                         graph::EdgeIsoMask edgemask)
        : cherrypicker(cp), small(f.GetGraph()), large(l), vmasks_large(vl),
          emasks_large(el), pmol(p), frag(f), has_mapping(false) {
      for (CMGV v : small.GetVertices())
        vmasks_small.emplace(v, v.GetIsomorphismMask() & vertmask);
      for (CMGE e : small.GetEdges())
        emasks_small.emplace(e, e.GetIsomorphismMask() & edgemask);
    }

    bool operator()(const CorrespondenceMap &map) override {
      using ConSym = graph::CMGVertex::ContractedSymmetry;
      has_mapping = true;
      std::vector<graph::MGVertex> frag_v, target_v;
      graph::CondensedMolecularGraph G = frag.GetGraph();
      graph::MolecularGraph molG = G.GetSuperGraph().GetMolecularGraph();
      Molecule fragMol = molG.GetMolecule();
      frag_v.reserve(molG.NumVertices());
      target_v.reserve(molG.NumVertices());
      std::vector<std::pair<size_t, size_t>> regions;
      for (auto &frag2target : map) {
        frag_v.emplace_back(frag2target.first.GetSource());
        target_v.emplace_back(frag2target.second.GetSource());
        size_t begin_size = frag_v.size();
        ConSym currentSym = ConSym::Hydrogen;
        for (auto &cv : frag2target.first.GetCondensedVertices()) {
          if (cv.first != currentSym) {
            regions.emplace_back(begin_size, frag_v.size());
            currentSym = cv.first;
            begin_size = frag_v.size();
          }
          frag_v.emplace_back(cv.second);
        }
        for (auto &cv : frag2target.second.GetContractedVertices())
          target_v.emplace_back(cv);
        regions.emplace_back(begin_size, frag_v.size());
      }

      RegionalPermutation permutation(frag_v.begin(), frag_v.end());
      for (std::pair<size_t, size_t> be : regions) {
        if (be.second - be.first < 2) continue;
        permutation.AddRegion(frag_v.begin() + be.first,
                              frag_v.begin() + be.second);
      }

      while (permutation()) {
        // Parameterise the atoms
        auto patms = frag.GetAtoms();
        for (size_t i = 0; i < frag_v.size(); ++i) {
          if (std::find(patms.begin(), patms.end(), frag_v[i]) == patms.end())
            continue;
          ParamAtom patm = pmol.GetAtom(target_v[i].GetAtom());
          patm.MappedWith(frag_v[i].GetAtom());
        }

        // Parameterise the bonds
        for (auto bnd : frag.GetBonds()) {
          graph::MGVertex v1 = bnd.first, v2 = bnd.second;
          if (!cherrypicker.GetBool(CPSet::AllowDanglingBonds) &&
              (std::find(patms.begin(), patms.end(), v1) == patms.end() ||
               std::find(patms.begin(), patms.end(), v2) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          Atom t1 = target_v[std::distance(frag_v.begin(), p1)].GetAtom();
          Atom t2 = target_v[std::distance(frag_v.begin(), p2)].GetAtom();
          ParamBond pbnd = pmol.GetBond(t1, t2);
          pbnd.MappedWith(fragMol.GetBond(v1.GetAtom(), v2.GetAtom()));
        }

        // Parameterise the angles
        for (auto ang : frag.GetAngles()) {
          graph::MGVertex v1 = ang.first, v2 = ang.second, v3 = ang.third;
          if (!cherrypicker.GetBool(CPSet::AllowDanglingAngles) &&
              (std::find(patms.begin(), patms.end(), v1) == patms.end() ||
               std::find(patms.begin(), patms.end(), v2) == patms.end() ||
               std::find(patms.begin(), patms.end(), v3) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          auto p3 = std::find(frag_v.begin(), frag_v.end(), v3);
          Atom t1 = target_v[std::distance(frag_v.begin(), p1)].GetAtom();
          Atom t2 = target_v[std::distance(frag_v.begin(), p2)].GetAtom();
          Atom t3 = target_v[std::distance(frag_v.begin(), p3)].GetAtom();
          ParamAngle pang = pmol.GetAngle(t1, t2, t3);
          pang.MappedWith(
              fragMol.GetAngle(v1.GetAtom(), v2.GetAtom(), v3.GetAtom()));
        }

        // Parameterise the dihedrals
        for (auto dhd : frag.GetDihedrals()) {
          graph::MGVertex v1 = dhd.first, v2 = dhd.second, v3 = dhd.third,
                          v4 = dhd.fourth;
          if (!cherrypicker.GetBool(CPSet::AllowDanglingDihedrals) &&
              (std::find(patms.begin(), patms.end(), v1) == patms.end() ||
               std::find(patms.begin(), patms.end(), v2) == patms.end() ||
               std::find(patms.begin(), patms.end(), v3) == patms.end() ||
               std::find(patms.begin(), patms.end(), v4) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          auto p3 = std::find(frag_v.begin(), frag_v.end(), v3);
          auto p4 = std::find(frag_v.begin(), frag_v.end(), v4);
          Atom t1 = target_v[std::distance(frag_v.begin(), p1)].GetAtom();
          Atom t2 = target_v[std::distance(frag_v.begin(), p2)].GetAtom();
          Atom t3 = target_v[std::distance(frag_v.begin(), p3)].GetAtom();
          Atom t4 = target_v[std::distance(frag_v.begin(), p4)].GetAtom();
          ParamDihedral pdhd = pmol.GetDihedral(t1, t2, t3, t4);
          pdhd.MappedWith(fragMol.GetDihedral(v1.GetAtom(), v2.GetAtom(),
                                              v3.GetAtom(), v4.GetAtom()));
        }
        if (!cherrypicker.GetBool(CPSet::ParameteriseFromAllPermutations))
          break;
      }
      return true;
    }

    bool operator()(const CMGV &vs, const CMGV &vl) override {
      return vmasks_small.at(vs) == vmasks_large.at(vl);
    }

    bool operator()(const CMGE &es, const CMGE &el) override {
      return emasks_small.at(es) == emasks_large.at(el);
    }
  };

  struct RICherryPickerMatcher : rilib::MatchListener {
    CherryPickerCallback &cb;
    RICherryPickerMatcher(CherryPickerCallback &_cb) : cb(_cb) {}
    virtual void match(int n, int *small_ids, int *large_ids) {
      CherryPickerCallback::CorrespondenceMap cp_map;
      for (int i = 0; i < n; ++i) {
        cp_map.emplace(cb.small.GetVertices()[small_ids[i]],
                       cb.large.GetVertices()[large_ids[i]]);
      }
      cb(cp_map);
    }
  };

  ParamMolecule CherryPicker::ParameteriseMolecule(Molecule &mol) {
    if (_libs.empty())
      throw std::runtime_error("No Athenaeums to parameterise from");
    graph::MolecularGraph G = mol.GetGraph();
    if (!G.IsConnected())
      throw std::runtime_error("CherryPicker requires a connected molecule");
    mol.PerceiveAngles();
    mol.PerceiveDihedrals();
    graph::CondensedMolecularGraph CMG = graph::Condense(G);
    ParamMolecule pmol(mol);

    // Populate the masks
    graph::VertexIsoMask vertmask;
    vertmask.reset();
    if (GetBool(CPSet::VertexElement)) vertmask |= graph::VertexIsoMask(0x7F);
    if (GetBool(CPSet::VertexFormalCharge))
      vertmask |= graph::VertexIsoMask(0x780);
    if (GetBool(CPSet::VertexCondensed)) {
      graph::VertexIsoMask tmp;
      tmp.from_uint64(0x1C03FFF800);
      vertmask |= tmp;
    }
    if (GetBool(CPSet::VertexCyclic))
      vertmask |= graph::VertexIsoMask(0xC000000);
    if (GetBool(CPSet::VertexStereochemistry))
      vertmask |= graph::VertexIsoMask(0x30000000);
    if (GetBool(CPSet::VertexAromaticity))
      vertmask |= graph::VertexIsoMask(0x40000000);
    if (GetBool(CPSet::VertexDegree)) {
      graph::VertexIsoMask tmp;
      tmp.from_uint64(0x380000000);
      vertmask |= tmp;
    }

    graph::EdgeIsoMask edgemask;
    edgemask.reset();
    if (GetBool(CPSet::EdgeBondOrder)) edgemask |= graph::EdgeIsoMask(7);
    if (GetBool(CPSet::EdgeStereochemistry)) edgemask |= graph::EdgeIsoMask(24);
    if (GetBool(CPSet::EdgeCyclic)) edgemask |= graph::EdgeIsoMask(96);
    //    if (GetBool(CPSet::EdgeAromaticity)) edgemask |=
    //    graph::EdgeIsoMask(128);
    if (GetBool(CPSet::EdgeDegree)) edgemask |= graph::EdgeIsoMask(16128);

    CherryPickerCallback::VertMasks vmasks;
    for (CMGV v : CMG.GetVertices())
      vmasks.emplace(v, vertmask & v.GetIsomorphismMask());
    CherryPickerCallback::EdgeMasks emasks;
    for (CMGE e : CMG.GetEdges())
      emasks.emplace(e, edgemask & e.GetIsomorphismMask());

    std::unique_ptr<rilib::Graph> CMG_ri;
    if (GetBool(CPSet::UseRISubgraphMatching))
      CMG_ri = CMGToRIGraph(CMG, edgemask, vertmask);

    // Run the matching
    for (Athenaeum &lib : _libs) {
      for (auto &g_frag : lib.GetFragments()) {
        // Initially all fragments are to be searched
        boost::dynamic_bitset<> fragments(g_frag.second.size());
        fragments.set();
        for (size_t pos = 0; pos < fragments.size();
             pos = fragments.find_next(pos)) {
          Fragment frag = g_frag.second[pos];
          if ((int32_t)frag.Size() < GetInt(CPSet::MinimumFragmentSize))
            continue;
          if (GetInt(CPSet::MaximumFragmentSize) > 0 &&
              (int32_t)frag.Size() > GetInt(CPSet::MaximumFragmentSize))
            continue;
          if (frag.GetGraph().NumVertices() > CMG.NumVertices()) continue;
          CherryPickerCallback callback(*this, CMG, vmasks, emasks, pmol, frag,
                                        vertmask, edgemask);
          graph::CondensedMolecularGraph FG = frag.GetGraph();

          if (!GetBool(CPSet::UseRISubgraphMatching)) {
            SubgraphIsomorphisms(FG, CMG, callback);
          } else {
            std::unique_ptr<rilib::Graph> FG_ri =
                CMGToRIGraph(FG, edgemask, vertmask);
            rilib::AttributeComparator *vert_compare =
                new Uint64AttrComparator();
            rilib::AttributeComparator *edge_compare =
                new Uint32AttrComparator();
            rilib::MatchListener *listener =
                new RICherryPickerMatcher(callback);
            rilib::MaMaConstrFirst *mama = new rilib::MaMaConstrFirst(*FG_ri);
            mama->build(*FG_ri);
            long tmp_1, tmp_2, tmp_3;
            // run the matching
            rilib::match(*CMG_ri, *FG_ri, *mama, *listener,
                         rilib::MATCH_TYPE::MT_INDSUB, *vert_compare,
                         *edge_compare, &tmp_1, &tmp_2, &tmp_3);

            delete mama;
            delete listener;
            delete vert_compare;
            delete edge_compare;
          }
          if (!callback.has_mapping) { fragments -= frag.GetSupersets(); }
        }
      }
      using ATSet = Athenaeum::Settings;
      pmol.ApplyParameteristion(lib.GetBool(ATSet::SelfConsistent));
    }
    return pmol;
  }

} // namespace indigox::algorithm
