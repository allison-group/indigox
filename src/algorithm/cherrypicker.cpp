#include <algorithm>
#include <map>
#include <vector>

#include <indigox/algorithm/cherrypicker.hpp>
#include <indigox/algorithm/graph/isomorphism.hpp>
#include <indigox/classes/athenaeum.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>
#include <indigox/graph/condensed.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/utils/combinatronics.hpp>

namespace indigox::algorithm {
  
  bool IXCherryPicker::Settings::AllowDanglingBonds = true;
  bool IXCherryPicker::Settings::AllowDanglingAngles = true;
  bool IXCherryPicker::Settings::AllowDanglingDihedrals = true;
  
  bool IXCherryPicker::AddAthenaeum(Athenaeum library) {
    if (library->GetForcefield() != _ff) return false;
    _libs.push_back(library);
    return true;
  }
  
  bool IXCherryPicker::RemoveAthenaeum(Athenaeum library) {
    auto pos = std::find(_libs.begin(), _libs.end(), library);
    if (pos != _libs.end()) _libs.erase(pos);
    return pos != _libs.end();
  }
  
  struct CherryPickerMappingCallback
  : public MappingCallback<graph::IXCondensedMolecularGraph> {
    using BaseType = MappingCallback<graph::IXCondensedMolecularGraph>;
    using CorrespondenceMap = BaseType::CorrespondenceMap;
    
    CherryPickerMappingCallback(ParamMolecule p, Fragment f)
    : pmol(p), frag(f) { }
    
    bool operator()(const CorrespondenceMap& map) override {
      using Settings = IXCherryPicker::Settings;
      using ConSym = graph::IXCMGVertex::ContractedSymmetry;
      std::vector<graph::MGVertex> frag_v, target_v;
      graph::CondensedMolecularGraph G = frag->GetGraph();
      graph::MolecularGraph molG = G->GetSource();
      Molecule fragMol = frag->GetMolecule();
      frag_v.reserve(molG->NumVertices());
      target_v.reserve(molG->NumVertices());
      std::vector<std::pair<size_t, size_t>> regions;
      for (auto& frag2target : map) {
        frag_v.emplace_back(frag2target.first->GetSource());
        target_v.emplace_back(frag2target.second->GetSource());
        size_t begin_size = frag_v.size();
        ConSym currentSym = ConSym::Hydrogen;
        for (auto& cv : frag2target.first->GetContractedVertices()) {
          if (cv.first != currentSym) {
            regions.emplace_back(begin_size, frag_v.size());
            currentSym = cv.first;
            begin_size = frag_v.size();
          }
          frag_v.emplace_back(cv.second);
        }
        for (auto& cv : frag2target.second->GetContractedVertices())
          target_v.emplace_back(cv.second);
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
        auto patms = frag->GetAtomVertices();
        for (size_ i = 0; i < frag_v.size(); ++i) {
          if (std::find(patms.begin(), patms.end(), frag_v[i]) == patms.end())
            continue;
          ParamAtom patm = pmol->GetAtom(target_v[i]->GetAtom());
          patm->MappedWith(frag_v[i]->GetAtom());
        }
        
        // Parameterise the bonds
        for (auto bnd : frag->GetBondVertices()) {
          graph::MGVertex v1 = bnd.first, v2 = bnd.second;
          if (!Settings::AllowDanglingBonds
              && (std::find(patms.begin(), patms.end(), v1) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v2) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          Atom t1 = target_v[std::distance(frag_v.begin(), p1)]->GetAtom();
          Atom t2 = target_v[std::distance(frag_v.begin(), p2)]->GetAtom();
          ParamBond pbnd = pmol->GetBond({t1,t2});
          pbnd->MappedWith(fragMol->GetBond(v1->GetAtom(), v2->GetAtom()));
        }
        
        // Parameterise the angles
        for (auto ang : frag->GetAngleVertices()) {
          graph::MGVertex v1 = ang.first, v2 = ang.second, v3 = ang.third;
          if (!Settings::AllowDanglingAngles
              && (std::find(patms.begin(), patms.end(), v1) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v2) == patms.end()
                  ||std::find(patms.begin(), patms.end(), v3) == patms.end()))
            break;
          auto p1 = std::find(frag_v.begin(), frag_v.end(), v1);
          auto p2 = std::find(frag_v.begin(), frag_v.end(), v2);
          auto p3 = std::find(frag_v.begin(), frag_v.end(), v3);
          Atom t1 = target_v[std::distance(frag_v.begin(), p1)]->GetAtom();
          Atom t2 = target_v[std::distance(frag_v.begin(), p2)]->GetAtom();
          Atom t3 = target_v[std::distance(frag_v.begin(), p3)]->GetAtom();
          ParamAngle pang = pmol->GetAngle({t1,t2,t3});
          pang->MappedWith(fragMol->GetAngle(v1->GetAtom(), v2->GetAtom(),
                                             v3->GetAtom()));
        }
        
        // Parameterise the dihedrals
        for (auto dhd : frag->GetDihedralVertices()) {
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
          Atom t1 = target_v[std::distance(frag_v.begin(), p1)]->GetAtom();
          Atom t2 = target_v[std::distance(frag_v.begin(), p2)]->GetAtom();
          Atom t3 = target_v[std::distance(frag_v.begin(), p3)]->GetAtom();
          Atom t4 = target_v[std::distance(frag_v.begin(), p4)]->GetAtom();
          ParamDihedral pdhd = pmol->GetDihedral({t1,t2,t3,t4});
          pdhd->MappedWith(fragMol->GetDihedral(v1->GetAtom(), v2->GetAtom(),
                                                v3->GetAtom(), v4->GetAtom()));
        }
      }
      return true;
    }
    
    ParamMolecule pmol;
    Fragment frag;
  };
  
  ParamMolecule IXCherryPicker::ParameteriseMolecule(Molecule mol) const {
    if (_libs.empty()) throw std::runtime_error("No Athenaeums to parameterise from");
    graph::CondensedMolecularGraph CMG = graph::CondenseMolecularGraph(mol->GetGraph());
    ParamMolecule pmol = std::make_shared<IXParamMolecule>(mol);
    for (Athenaeum lib : _libs) {
      for (auto& g_frag : lib->GetFragments()) {
        for (Fragment frag : g_frag.second) {
          CherryPickerMappingCallback callback(pmol, frag);
          SubgraphIsomorphisms(frag->GetGraph(), CMG, callback);
        }
      }
      pmol->ApplyParameteristion(lib->IsManualSelfConsistent());
    }
    return pmol;
  }
  
}
