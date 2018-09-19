#include <boost/math/special_functions/relative_difference.hpp>

#include <indigox/classes/angle.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/classes/dihedral.hpp>
#include <indigox/classes/forcefield.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/parameterised.hpp>

namespace indigox {
  void IXParamAtom::MappedWith(Atom mapped) {
    if (_applied) return;
    FFAtom t = mapped->GetType();
    auto t_pos = _counts.find(t);
    if (t_pos == _counts.end()) _counts.emplace(t, 1);
    else ++t_pos->second;
    _charges.push_back(mapped->GetPartialCharge());
    _atms.emplace_back(mapped);
  }
  
  void IXParamAtom::ApplyParameterisation(bool self_consistent) {
    if (_atms.empty()) return;  // Only parametrise when parameters to apply
    if (_applied) return;  // Only parameterise if not done so already
    float_ mean = MeanCharge();
    if (self_consistent) {
      if (_counts.size() > 1)
        throw std::runtime_error("Types not self-consistent");
      if (boost::math::relative_difference(mean, MeadianCharge()) > 1e-10)
        throw std::runtime_error("Charges mean/median not equal");
      if (boost::math::relative_difference(0.0, StandardDeviationCharge()) > 1e-10)
        throw std::runtime_error("Charge stddev not 0");
    }
    Atom atm = _atm.lock();
    if (!atm) throw std::runtime_error("Mapped atom missing");
    atm->SetType(GetMostCommonType());
    atm->SetPartialCharge(mean);
    _applied = true;
  }
  
  void IXParamBond::MappedWith(Bond mapped) {
    if (_applied) return;
    FFBond t = mapped->GetType();
    FFBond t2 = t->GetLinkedType();
    auto t_pos = _counts.find(t);
    auto t2_pos = _counts.find(t2);
    auto end = _counts.end();
    
    if (t_pos == end && t2_pos == end) _counts.emplace(t, 1);
    else if (t2 && t2_pos != end) ++t2_pos->second;
    else ++t_pos->second;
    
    _bnds.emplace_back(mapped);
  }
  
  void IXParamBond::ApplyParameterisation(bool self_consistent) {
    if (_bnds.empty()) return;  // Only parametrise when parameters to apply
    if (_applied) return;  // Only parameterise if not done so already
    if (self_consistent && _counts.size() > 1)
      throw std::runtime_error("Types not self-consistent");
    Bond bnd = _bnd.lock();
    if (!bnd) throw std::runtime_error("Mapped bond missing");
    bnd->SetType(GetMostCommonType());
    _applied = true;
  }
  
  void IXParamAngle::MappedWith(Angle mapped) {
    if (_applied) return;
    FFAngle t = mapped->GetType();
    FFAngle t2 = t->GetLinkedType();
    auto t_pos = _counts.find(t);
    auto t2_pos = _counts.find(t2);
    auto end = _counts.end();
    
    if (t_pos == end && t2_pos == end) _counts.emplace(t, 1);
    else if (t2 && t2_pos != end) ++t2_pos->second;
    else ++t_pos->second;
    
    _angs.emplace_back(mapped);
  }
  
  void IXParamAngle::ApplyParameterisation(bool self_consistent) {
    if (_angs.empty()) return;
    if (_applied) return;
    if (self_consistent && _counts.size() > 1)
      throw std::runtime_error("Types not self-consistent");
    Angle ang = _ang.lock();
    if (!ang) throw std::runtime_error("Mapped angle missing");
    ang->SetType(GetMostCommonType());
    _applied = true;
  }
  
  void IXParamDihedral::MappedWith(Dihedral mapped) {
    if (_applied) return;
    FFDihedral t = mapped->GetType();
    auto t_pos = _counts.find(t);
    if (t_pos == _counts.end()) _counts.emplace(t, 1);
    else ++t_pos->second;
    _dhds.emplace_back(mapped);
  }
  
  void IXParamDihedral::ApplyParameterisation(bool self_consistent) {
    if (_dhds.empty()) return;
    if (_applied) return;
    if (self_consistent && _counts.size() > 1)
      throw std::runtime_error("Types not self-consistent");
    Dihedral dhd = _dhd.lock();
    if (!dhd) throw std::runtime_error("Mapped dihedral missing");
    dhd->SetType(GetMostCommonType());
    _applied = true;
  }
  
  IXParamMolecule::IXParamMolecule(Molecule mol) : _mol(mol) {
    for (auto it = mol->GetAtoms(); it.first != it.second; ++it.first)
      _atms.emplace(*it.first, std::make_shared<IXParamAtom>(*it.first));
    for (auto it = mol->GetBonds(); it.first != it.second; ++it.first) {
      PBond atms = (*it.first)->GetAtoms();
      _bnds.emplace(atms, std::make_shared<IXParamBond>(atms, (*it.first)));
    }
    for (auto it = mol->GetAngles(); it.first != it.second; ++it.first) {
      PAngle atms = (*it.first)->GetAtoms();
      _angs.emplace(atms, std::make_shared<IXParamAngle>(atms, (*it.first)));
    }
    for (auto it = mol->GetDihedrals(); it.first != it.second; ++it.first) {
      PDihedral atms = (*it.first)->GetAtoms();
      _dhds.emplace(atms, std::make_shared<IXParamDihedral>(atms, (*it.first)));
    }
  }
  
  ParamAtom IXParamMolecule::GetAtom(Atom atm) const {
    auto pos = _atms.find(atm);
    return pos == _atms.end() ? ParamAtom() : pos->second;
  }
  
  ParamBond IXParamMolecule::GetBond(Bond bnd) const {
    return GetBond(bnd->GetAtoms());
  }
  
  ParamBond IXParamMolecule::GetBond(IXParamMolecule::PBond atms) const {
    auto pos = _bnds.find(atms);
    if (pos != _bnds.end()) return pos->second;
    pos = _bnds.find(std::make_pair(atms.second, atms.first));
    return pos == _bnds.end() ? ParamBond() : pos->second;
  }
  
  ParamAngle IXParamMolecule::GetAngle(Angle ang) const {
    return GetAngle(ang->GetAtoms());
  }
  
  ParamAngle IXParamMolecule::GetAngle(IXParamMolecule::PAngle atms) const {
    auto pos = _angs.find(atms);
    if (pos != _angs.end()) return pos->second;
    pos = _angs.find(stdx::make_triple(atms.third, atms.second, atms.first));
    return pos == _angs.end() ? ParamAngle() : pos->second;
  }
  
  ParamDihedral IXParamMolecule::GetDihedral(Dihedral dhd) {
    return GetDihedral(dhd->GetAtoms());
  }
  
  ParamDihedral IXParamMolecule::GetDihedral(IXParamMolecule::PDihedral atms) {
    auto pos = _dhds.find(atms);
    if (pos != _dhds.end()) return pos->second;
    pos = _dhds.find(stdx::make_quad(atms.fourth, atms.third, atms.second, atms.first));
    if (pos != _dhds.end()) return pos->second;
    //! \todo Create a new dihedral in the molecule
    ParamDihedral dhd = std::make_shared<IXParamDihedral>(atms, Dihedral());
    _dhds.emplace(atms, dhd);
    return dhd;
  }
  
  void IXParamMolecule::ApplyParameteristion(bool self_consistent) {
    for (auto& atm : _atms) atm.second->ApplyParameterisation(self_consistent);
    for (auto& bnd : _bnds) bnd.second->ApplyParameterisation(self_consistent);
    for (auto& ang : _angs) ang.second->ApplyParameterisation(self_consistent);
    for (auto& dhd : _dhds) dhd.second->ApplyParameterisation(self_consistent);
  }
}
