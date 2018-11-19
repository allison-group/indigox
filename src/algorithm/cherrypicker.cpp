#include <algorithm>
#include <map>
#include <vector>

#include <indigox/algorithm/cherrypicker.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/combinatronics.hpp>
#include <indigox/algorithm/graph/isomorphism.hpp>

namespace indigox::algorithm {
  using VParam = CherryPicker::VertexParameters;
  using EParam = CherryPicker::EdgeParameters;
  using CPSet = CherryPicker::Settings;
  
  bool CPSet::AllowDanglingBonds = true;
  bool CPSet::AllowDanglingAngles = true;
  bool CPSet::AllowDanglingDihedrals = true;
  VParam CPSet::VertexMapping = (VParam::ElementType | VParam::FormalCharge
                                 | VParam::CondensedVertices | VParam::Degree);
  EParam CPSet::EdgeMapping = EParam::BondOrder | EParam::Degree;
  uint32_t CPSet::MinimumFragmentSize = 4;
  
  CherryPicker::CherryPicker(const Forcefield& ff) : _ff(ff) { }
  
  bool CherryPicker::AddAthenaeum(const Athenaeum& library) {
    if (library.GetForcefield() != _ff) return false;
    _libs.push_back(library);
    return true;
  }
  
  bool CherryPicker::RemoveAthenaeum(const Athenaeum& library) {
    auto pos = std::find(_libs.begin(), _libs.end(), library);
    if (pos != _libs.end()) _libs.erase(pos);
    return pos != _libs.end();
  }
  
  using CMGV = graph::CMGVertex;
  using CMGE = graph::CMGEdge;
  using CMGS = graph::sCondensedMolecularGraph;
  using U = graph::Undirected;
  using GL = graph::GraphLabel;
  
  struct CherryPickerCallback
  : public CMGCallback {
    using GraphType = graph::CondensedMolecularGraph;
    using BaseType = CMGCallback;
    using CorrespondenceMap = BaseType::CorrespondenceMap;
    
    using VertMasks = eastl::vector_map<CMGV, graph::VertexIsoMask>;
    using EdgeMasks = eastl::vector_map<CMGE, graph::EdgeIsoMask>;
    
    GraphType& small;
    GraphType& large;
    VertMasks& vmasks_large;
    EdgeMasks& emasks_large;
    VertMasks vmasks_small;
    EdgeMasks emasks_small;
    ParamMolecule& pmol;
    Fragment& frag;
    
    CherryPickerCallback(GraphType& l, VertMasks& vl, EdgeMasks& el,
                         ParamMolecule& p, Fragment& f,
                         graph::VertexIsoMask vertmask, graph::EdgeIsoMask edgemask)
    : small(f.GetGraph()), large(l), vmasks_large(vl), emasks_large(el), pmol(p),
    frag(f) {
      for (CMGV v : small.GetVertices())
        vmasks_small.emplace(v, v.GetIsomorphismMask() & vertmask);
      for (CMGE e : small.GetEdges())
        emasks_small.emplace(e, e.GetIsomorphismMask() & edgemask);
      
    }
    
    bool operator()(const CorrespondenceMap& map) override {
      using Settings = CherryPicker::Settings;
      using ConSym = graph::CMGVertex::ContractedSymmetry;
      std::vector<graph::MGVertex> frag_v, target_v;
      graph::CondensedMolecularGraph& G = frag.GetGraph();
      graph::MolecularGraph& molG = G.GetSuperGraph().GetMolecularGraph();
      Molecule& fragMol = molG.GetMolecule();
      frag_v.reserve(molG.NumVertices());
      target_v.reserve(molG.NumVertices());
      std::vector<std::pair<size_t, size_t>> regions;
      for (auto& frag2target : map) {
        frag_v.emplace_back(frag2target.first.GetSource());
        target_v.emplace_back(frag2target.second.GetSource());
        size_t begin_size = frag_v.size();
        ConSym currentSym = ConSym::Hydrogen;
        for (auto& cv : frag2target.first.GetCondensedVertices()) {
          if (cv.first != currentSym) {
            regions.emplace_back(begin_size, frag_v.size());
            currentSym = cv.first;
            begin_size = frag_v.size();
          }
          frag_v.emplace_back(cv.second);
        }
        for (auto& cv : frag2target.second.GetContractedVertices())
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
          if (!Settings::AllowDanglingBonds
              && (std::find(patms.begin(), patms.end(), v1) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v2) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          sAtom t1 = target_v[std::distance(frag_v.begin(), p1)].GetAtom().shared_from_this();
          sAtom t2 = target_v[std::distance(frag_v.begin(), p2)].GetAtom().shared_from_this();
          ParamBond pbnd = pmol.GetBond({t1,t2});
          pbnd.MappedWith(fragMol.GetBond(v1.GetAtom(), v2.GetAtom()));
        }
        
        // Parameterise the angles
        for (auto ang : frag.GetAngles()) {
          graph::MGVertex v1 = ang.first, v2 = ang.second, v3 = ang.third;
          if (!Settings::AllowDanglingAngles
              && (std::find(patms.begin(), patms.end(), v1) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v2) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v3) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          auto p3 = std::find(frag_v.begin(), frag_v.end(), v3);
          sAtom t1 = target_v[std::distance(frag_v.begin(), p1)].GetAtom().shared_from_this();
          sAtom t2 = target_v[std::distance(frag_v.begin(), p2)].GetAtom().shared_from_this();
          sAtom t3 = target_v[std::distance(frag_v.begin(), p3)].GetAtom().shared_from_this();
          ParamAngle pang = pmol.GetAngle({t1,t2,t3});
          pang.MappedWith(fragMol.GetAngle(v1.GetAtom(), v2.GetAtom(), v3.GetAtom()));
        }
        
        // Parameterise the dihedrals
        for (auto dhd : frag.GetDihedrals()) {
          graph::MGVertex v1 = dhd.first, v2 = dhd.second, v3 = dhd.third, v4 = dhd.fourth;
          if (!Settings::AllowDanglingDihedrals
              && (std::find(patms.begin(), patms.end(), v1) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v2) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v3) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v4) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          auto p3 = std::find(frag_v.begin(), frag_v.end(), v3);
          auto p4 = std::find(frag_v.begin(), frag_v.end(), v4);
          sAtom t1 = target_v[std::distance(frag_v.begin(), p1)].GetAtom().shared_from_this();
          sAtom t2 = target_v[std::distance(frag_v.begin(), p2)].GetAtom().shared_from_this();
          sAtom t3 = target_v[std::distance(frag_v.begin(), p3)].GetAtom().shared_from_this();
          sAtom t4 = target_v[std::distance(frag_v.begin(), p4)].GetAtom().shared_from_this();
          ParamDihedral pdhd = pmol.GetDihedral({t1,t2,t3,t4});
          pdhd.MappedWith(fragMol.GetDihedral(v1.GetAtom(), v2.GetAtom(),
                                              v3.GetAtom(), v4.GetAtom()));
        }
      }
      return true;
    }
    
    bool operator()(const CMGV& vs, const CMGV& vl) override {
      return vmasks_small.at(vs) == vmasks_large.at(vl);
    }
    
    bool operator()(const CMGE& es, const CMGE& el) override {
      return emasks_small.at(es) == emasks_large.at(el);
    }
    
  };
  
  ParamMolecule CherryPicker::ParameteriseMolecule(Molecule& mol) {
    if (_libs.empty())
      throw std::runtime_error("No Athenaeums to parameterise from");
    if (!mol.GetGraph().IsConnected())
      throw std::runtime_error("CherryPicker requires a connected molecule");
    graph::sCondensedMolecularGraph CMG = graph::Condense(mol.GetGraph());
    ParamMolecule pmol(mol);
    
    // Populate the masks
    graph::VertexIsoMask vertmask; vertmask.reset();
    if ((CPSet::VertexMapping & VParam::ElementType) != VParam::None)
      vertmask |= graph::VertexIsoMask(0x7F);
    if ((CPSet::VertexMapping & VParam::FormalCharge) != VParam::None)
      vertmask |= graph::VertexIsoMask(0x780);
    if ((CPSet::VertexMapping & VParam::CondensedVertices) != VParam::None)
      vertmask |= graph::VertexIsoMask(0x3FFF800);
    if ((CPSet::VertexMapping & VParam::CyclicNature) != VParam::None)
      vertmask |= graph::VertexIsoMask(0xC000000);
    if ((CPSet::VertexMapping & VParam::Stereochemistry) != VParam::None)
      vertmask |= graph::VertexIsoMask(0x30000000);
    if ((CPSet::VertexMapping & VParam::Aromaticity) != VParam::None)
      vertmask |= graph::VertexIsoMask(0x40000000);
    if ((CPSet::VertexMapping & VParam::Degree) != VParam::None) {
      graph::VertexIsoMask tmp; tmp.from_uint64(0x380000000);
      vertmask |= tmp;
    }
    
    graph::EdgeIsoMask edgemask; edgemask.reset();
    if ((CPSet::EdgeMapping & EParam::BondOrder) != EParam::None)
      edgemask |= graph::EdgeIsoMask(7);
    if ((CPSet::EdgeMapping & EParam::Stereochemistry) != EParam::None)
      edgemask |= graph::EdgeIsoMask(24);
    if ((CPSet::EdgeMapping & EParam::CyclicNature) != EParam::None)
      edgemask |= graph::EdgeIsoMask(96);
    if ((CPSet::EdgeMapping & EParam::Aromaticity) != EParam::None)
      edgemask |= graph::EdgeIsoMask(128);
    if ((CPSet::EdgeMapping & EParam::Degree) != EParam::None)
      edgemask |= graph::EdgeIsoMask(16128);
    
    CherryPickerCallback::VertMasks vmasks;
    for (CMGV v : CMG->GetVertices())
      vmasks.emplace(v, vertmask & v.GetIsomorphismMask());
    CherryPickerCallback::EdgeMasks emasks;
    for (CMGE e : CMG->GetEdges())
      emasks.emplace(e, edgemask & e.GetIsomorphismMask());
    
    // Run the matching
    for (Athenaeum& lib : _libs) {
      for (auto& g_frag : lib.GetFragments()) {
        for (Fragment frag : g_frag.second) {
          if (frag.Size() < CPSet::MinimumFragmentSize) continue;
          if (frag.GetGraph().NumVertices() > CMG->NumVertices()) continue;
          CherryPickerCallback callback(*CMG, vmasks, emasks, pmol, frag,
                                        vertmask, edgemask);
          SubgraphIsomorphisms(frag.GetGraph(), *CMG, callback);
        }
      }
      pmol.ApplyParameteristion(lib.IsSelfConsistent());
    }
    return pmol;
  }
  
}
