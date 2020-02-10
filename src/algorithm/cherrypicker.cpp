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
#include <indigo-bondorder/indigo-bondorder.hpp>

#include <algorithm>
#include <map>
#include <vector>

namespace indigox::algorithm {

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
    SetBool(CPSet::CalculateElectrons);

    SetInt(CPSet::MinimumFragmentSize, 4);
    SetInt(CPSet::MaximumFragmentSize, -1);
    SetInt(CPSet::ElectronMethod, 2); // Default to FPT for high accuracy and speed
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
    std::cout << "Adding new Athenaeum..." << std::endl;
    if (library.GetForcefield() != _ff) return false;
    _libs.push_back(library);
    return true;
  }

  bool CherryPicker::RemoveAthenaeum(Athenaeum &library) {
    std::cout << "Removing Athenaeum..." << std::endl;
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
    std::cout << "Parameterising molecule " << mol.GetName() << "." << std::endl;

    if (_libs.empty())
      throw std::runtime_error("No Athenaeums to parameterise from");
    graph::MolecularGraph G = mol.GetGraph();
    if (!G.IsConnected())
      throw std::runtime_error("CherryPicker requires a connected molecule");
    mol.PerceiveAngles();
    mol.PerceiveDihedrals();

    if (GetBool(CPSet::CalculateElectrons)) {
      mol.PerceiveElectrons(GetInt(CPSet::ElectronMethod), GetBool(CPSet::NoInput));
    }

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
    if (GetInt(CPSet::ChargeRounding) > 3) {
      ParamMolecule::charge_rounding = 1;
      for (int32_t i = 0; i < GetInt(CPSet::ChargeRounding); ++i)
        ParamMolecule::charge_rounding *= 10;
    } else {
      ParamMolecule::charge_rounding = 1000;
    }
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
            std::unique_ptr<rilib::Graph> FG_ri = CMGToRIGraph(FG, edgemask, vertmask);
            std::unique_ptr<rilib::AttributeComparator> vert_compare = std::make_unique<Uint64AttrComparator>();
            std::unique_ptr<rilib::AttributeComparator> edge_compare = std::make_unique<Uint32AttrComparator>();
            std::unique_ptr<rilib::MatchListener> listener = std::make_unique<RICherryPickerMatcher>(callback);
            std::unique_ptr<rilib::MaMaConstrFirst> mama = std::make_unique<rilib::MaMaConstrFirst>(*FG_ri);
            mama->build(*FG_ri);
            long tmp_1, tmp_2, tmp_3;
            // run the matching
            rilib::match(*CMG_ri, *FG_ri, *mama, *listener,
                         rilib::MATCH_TYPE::MT_INDSUB, *vert_compare,
                         *edge_compare, &tmp_1, &tmp_2, &tmp_3);
          }
          if (!callback.has_mapping) { fragments -= frag.GetSupersets(); }
        }
      }
      using ATSet = Athenaeum::Settings;
      pmol.ApplyParameteristion(lib.GetBool(ATSet::SelfConsistent));
    }
    
    // Redistribute any excess charge, but only if all atoms have been mapped
    bool redistribute = true;
    for (ParamAtom patm : pmol.GetAtoms()) {
      if (patm.GetMappedCharges().empty()) {
        redistribute = false;
        break;
      }
    }
    
    if (redistribute) {
      std::vector<ParamAtom> addable_atoms = pmol.GetChargeAddableAtoms();
      std::sort(addable_atoms.begin(), addable_atoms.end(), [](ParamAtom& a, ParamAtom& b) { return a.MeanCharge() < b.MeanCharge(); });
      
      
      double target_charge = mol.GetMolecularCharge();
      double total_charge = 0.;
      for (Atom atm : mol.GetAtoms()) total_charge += atm.GetPartialCharge();
      double to_add = target_charge - total_charge;
      uint64_t count = (uint64_t)abs(round(to_add * ParamMolecule::charge_rounding));
      if (count && addable_atoms.empty()) {
        std::cout << "WARNING: Total charge does not match target charge but no atoms are available for charge redistribution.\n";
        return pmol;
      }
      if (abs(to_add) > 0.1) std::cout << "WARNING: CherryPicker redistributing a large charge imbalance: " << to_add << "\n";
      
      if (to_add < 0) {
        double charge_delta = -1. / ParamMolecule::charge_rounding;
        for (int64_t pos = 0; count; count -= 1, pos += 1) {
          if (pos == (int64_t)addable_atoms.size()) pos = 0;
          addable_atoms[pos].AddRedistributedCharge(charge_delta);
        }
      } else {
        double charge_delta = 1. / ParamMolecule::charge_rounding;
        for (int64_t pos = addable_atoms.size() - 1; count; count -= 1, pos -= 1) {
          if (!pos) pos = addable_atoms.size() - 1;
          addable_atoms[pos].AddRedistributedCharge(charge_delta);
        }
      }
    } else std::cout << "WARNING: Not all atoms mapped so charge cannot be redistributed.\n";

    std::cout << "Finished parameterising molecule " << mol.GetName() << ".\n" << std::endl;
    
    return pmol;
  }

} // namespace indigox::algorithm
