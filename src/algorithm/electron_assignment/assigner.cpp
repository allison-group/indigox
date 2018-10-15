#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <set>

#include <EASTL/vector_set.h>

#include <boost/algorithm/string.hpp>

#include <indigox/algorithm/electron_assignment/astar_optimisation.hpp>
#include <indigox/algorithm/electron_assignment/assigner.hpp>
#include <indigox/algorithm/electron_assignment/local_optimisation.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>
#include <indigox/classes/molecule.hpp>
#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/filereader.hpp>

namespace indigox::algorithm {
  using Option = utils::Option;
  using algo = AssignerAlgorithm;
  using settings = IXElectronAssigner::Settings;
  using AAlgo = IXElectronAssigner::AssignAlgorithm;
  
  // Default states for the ElectronAssigner settings
  AssignerAlgorithm settings::Algorithm = AssignerAlgorithm::AStar;
  Option settings::ElectronPairs = Option::Default;
  Option settings::ChargedCarbon = Option::Default;
  Option settings::Preassign = Option::Default;
  std::string settings::ScoreFile = "assignment_scores.dat";
  score_t settings::Infinity = std::numeric_limits<score_t>::max();
  uint8_t settings::MaxBondOrder = 3;
  int16_t settings::MaxChargeMagnitude = -1;
  uint16_t settings::MaxNumResults = 64;
  eastl::vector_set<Element> settings::AllowedElements = {
    GetPeriodicTable().GetElement("H"),
    GetPeriodicTable().GetElement("C"),
    GetPeriodicTable().GetElement("N"),
    GetPeriodicTable().GetElement("O"),
    GetPeriodicTable().GetElement("S"),
    GetPeriodicTable().GetElement("P"),
    GetPeriodicTable().GetElement("F"),
    GetPeriodicTable().GetElement("Cl"),
    GetPeriodicTable().GetElement("Br")
  };
  
  AAlgo::AssignAlgorithm(const ScoreTable& t)
  : _mol(), _g(), _locs(), _num_e(0), _opts(0), _inf(settings::Infinity),
  _min_score(_inf), _limit(_inf), _results(),
  _table(t), _max_charge(settings::MaxChargeMagnitude),
  _max_results(settings::MaxNumResults), _previous_mask(0)
  {
    if (settings::ElectronPairs == Option::Yes
        || settings::ElectronPairs == Option::Default)
      _opts.set(__pairs);
    if (settings::ChargedCarbon == Option::Yes
        || settings::ChargedCarbon == Option::Default)
      _opts.set(__charged_c);
    if (settings::Preassign == Option::Yes
        || settings::Preassign == Option::Default)
      _opts.set(__preplace);
  }
  
  void AAlgo::Initalise(const Molecule& mol) {
    _mol = mol;
    _g.reset();
    
    // Make the assignment graph
    _g = std::make_shared<graph::IXAssignmentGraph>(mol->GetGraph());
    if (_opts[__preplace]) _g->PreassignElectrons();
    
    // Calculate number of electrons to assign
    int count = -mol->GetMolecularCharge();
    for (auto its = mol->GetAtomIters(); its.first != its.second; ++its.first)
      count += (*its.first)->GetElement().GetValenceElectronCount();
    for (auto its = _g->GetVertices(); its.first != its.second; ++its.first)
      count -= (*its.first)->GetPreAssignedCount();
    
    // Set if using pairs or not
    if (settings::ElectronPairs == Option::Auto && count % 2)
      _opts.reset(__pairs);
    else if (settings::ElectronPairs == Option::Auto)
      _opts.set(__pairs);
    
    // Pairs sanity checl
    if (_opts[__pairs] && count % 2)
      throw std::runtime_error("Odd number of electrons when using pairs");
    else if (_opts[__pairs]) count /= 2;
    
    uint32_t delta = 2;
    if (!_opts[__pairs]) delta = 1;
    
    // Range sanity checks
    if (count < 0)
      throw std::runtime_error("Negative electrons to assign");
    _num_e = static_cast<size_t>(count);
    
    // Figure out where to place them Vertices
    indigox::graph::IXMolecularGraph& G = *mol->GetGraph();
    for (auto its = G.GetVertices(); its.first != its.second; ++its.first) {
      indigox::graph::MGVertex v = *its.first;
      int octet;
      if (G.Degree(v) > 2)
        octet = (v)->GetAtom()->GetElement().GetHypervalentOctet();
      else
        octet = (v)->GetAtom()->GetElement().GetOctet();
      
      int bonded = 2 * G.Degree(v);
      int missing = octet - bonded;
      
      indigox::graph::AGVertex av = _g->GetVertex(v);
      missing -= av->GetPreAssignedCount();
      
      // Missing sanity checks
      if (_opts[__pairs] && missing % 2)
        throw std::runtime_error("Odd number of electrons when using pairs");
      else if (missing < 0)
        throw std::runtime_error("Negative electrons missing");
      
      while (missing) {
        _locs.emplace_back(av);
        missing -= delta;
      }
    }
    
    // Figure out where to place them edges
    for (auto its = G.GetEdges(); its.first != its.second; ++its.first) {
      indigox::graph::MGEdge e = *its.first;
      indigox::graph::MGVertex u = G.GetSource(e), v = G.GetTarget(e);
      int u_oct;
      if (G.Degree(u) > 2)
        u_oct = u->GetAtom()->GetElement().GetHypervalentOctet();
      else
        u_oct = u->GetAtom()->GetElement().GetOctet();
      int u_bonded = 2 * G.Degree(u);
      int u_missing = u_oct - u_bonded;
      
      int v_oct;
      if (G.Degree(v) > 2)
        v_oct = v->GetAtom()->GetElement().GetHypervalentOctet();
      else
        v_oct = v->GetAtom()->GetElement().GetOctet();
      int v_bonded = 2 * G.Degree(v);
      int v_missing = v_oct - v_bonded;
      
      indigox::graph::AGVertex av = _g->GetVertex(e);
      uint32_t missing = 2 * settings::MaxBondOrder;
      if (missing < av->GetPreAssignedCount())
        throw std::runtime_error("Too many pre-assigned electrons");
      missing -= av->GetPreAssignedCount();
      while (missing > 0 && u_missing > 0 && v_missing > 0) {
        _locs.emplace_back(av);
        missing -= delta;
        u_missing -= delta;
        v_missing -= delta;
      }
    }
    // Final range sanity check
    if (_num_e > _locs.size())
      throw std::runtime_error("More to assign than places to put.");
    
    // Set the initial bitset value
    _previous_mask = AssignMask(_locs.size());
    // Set that the algorithm has been initalised.
    _opts.set(__initalised);
  }

  graph::AssignmentGraph AAlgo::GetAssignmentGraph() {
    if (_mol.expired()) return graph::AssignmentGraph();
    else if (!_g)
      _g = std::make_shared<graph::IXAssignmentGraph>(_mol.lock()->GetGraph());
    return _g;
  }
  
  AssignMask AAlgo::CalculateUpperLimit() {
    if (!IsInitalised())
      throw std::runtime_error("Requires initalised algorithm");
    AssignMask initial = _previous_mask;
    initial.reset();
    
    std::vector<size_t> all_location(_locs.size());
    std::iota(all_location.begin(), all_location.end(), 0);
    
    while (initial.count() < _num_e) {
      score_t min_nbr_score = _inf;
      size_t min_nbr_pos = all_location.front();
      for (size_t i : all_location) {
        initial.set(i);
        score_t current_score = CalculateAssignmentScore(initial);
        if (current_score < min_nbr_score)
          std::tie(min_nbr_score, min_nbr_pos) = {current_score, i};
        initial.reset(i);
      }
      initial.set(min_nbr_pos);
      all_location.erase(std::remove(all_location.begin(), all_location.end(),
                                     min_nbr_pos));
    }
    _limit = CalculateAssignmentScore(initial);
    if (_limit < _inf) ++_limit;
    return initial;
  }
  
  score_t AAlgo::CalculateVertexScore(const graph::AGVertex &v) const {
    if (!IsInitalised())
      throw std::runtime_error("Requires initalised algorithm");
    size_t occurances = std::count(_locs.begin(), _locs.end(), v);
    key_t k = 0;
    if (v->IsVertexMapped()) {
      const Element& e = v->GetSourceVertex()->GetAtom()->GetElement();
      char fc = e.GetValenceElectronCount() - v->GetTotalAssigned();
      uint8_t valence = v->GetTotalAssigned();
      size_t nbr_count = 0;
      auto nbrs_it = _g->GetNeighbours(v);
      for (; nbrs_it.first != nbrs_it.second; ++nbrs_it.first) {
        fc -= (*nbrs_it.first)->GetTotalAssigned() / 2;
        valence += (*nbrs_it.first)->GetTotalAssigned();
        nbr_count += std::count(_locs.begin(), _locs.end(), *nbrs_it.first);
      }
      
      // Allow charged carbons check
      if (!_opts[__charged_c] && e == "C" && fc != 0) return _inf;
      // Highest magnitude charge check
      if (_max_charge >= 0 && std::abs(fc) > _max_charge) return _inf;
      // See if we can set the score to 0. Can do so when both the number of
      // occurances of the vertex and all its neighbours is 0
      if (occurances == 0 && nbr_count == 0) return 0;
      // Valence check (valence must be less than or equal target octet
      if ((_g->Degree(v) > 2 && valence > e.GetHypervalentOctet())
          || (valence > e.GetOctet())) return _inf;
      
      // Create the key
      k = e.GetAtomicNumber() + (std::abs(fc) << 8);
      if (fc < 0) k += 1 << 15;
    } else {
      // Can set the score to 0 when no occurances of the vertex
      if (occurances == 0) return 0;
      
      auto nbrs_it = _g->GetNeighbours(v);
      graph::AGVertex u1 = *nbrs_it.first++;
      graph::AGVertex u2 = *nbrs_it.first;
      
      uint8_t valence_u1 = u1->GetTotalAssigned();
      uint8_t valence_u2 = u2->GetTotalAssigned();
      for (auto it = _g->GetNeighbours(u1); it.first != it.second; ++it.first)
        valence_u1 += (*it.first)->GetTotalAssigned();
      for (auto it = _g->GetNeighbours(u2); it.first != it.second; ++it.first)
        valence_u2 += (*it.first)->GetTotalAssigned();
      Element e1 = u1->GetSourceVertex()->GetAtom()->GetElement();
      Element e2 = u1->GetSourceVertex()->GetAtom()->GetElement();
      // Valence check of vertices making up edge vertex
      if ((_g->Degree(u1) > 2 && valence_u1 > e1.GetHypervalentOctet())
          || (valence_u1 > e1.GetOctet())) return _inf;
      if ((_g->Degree(u2) > 2 && valence_u2 > e2.GetHypervalentOctet())
          || (valence_u2 > e2.GetOctet())) return _inf;
      
      // Create the key
      k = e1.GetAtomicNumber() + (e2.GetAtomicNumber() << 8);
      k += v->GetTotalAssigned() << 20;
    }
    // Find the key score
    auto pos = _table.find(k);
    return pos == _table.end() ? _inf : pos->second;
  }
  
  score_t AAlgo::CalculateAssignmentScore(const AssignMask &a) {
    SetAssignment(a);
    
    score_t score = 0;
    for (auto vs = _g->GetVertices(); vs.first != vs.second; ++vs.first) {
      score_t vscore = CalculateVertexScore(*vs.first);
      if (vscore == _inf) return _inf;
      score += vscore;
    }
    return score;
  }
  
  void AAlgo::SetAssignment(const AssignMask &a) { 
    AssignMask changed = _previous_mask ^ a;
    AssignMask lost = _previous_mask & a;
    uint32_t delta = 2;
    if (!_opts[__pairs]) delta = 1;
    size_t i = changed.find_first();
    while (i < changed.size()) {
      graph::AGVertex v= _locs[i];
      if (lost[i]) v->SetAssignedCount(v->GetAssignedCount() - delta);
      else v->SetAssignedCount(v->GetAssignedCount() + delta);
      i = changed.find_next(i);
    }
    _previous_mask = a;
  }
  
  bool AAlgo::ApplyAssignment(size_t idx) {
    if (!IsInitalised())
      throw std::runtime_error("Requires initalistion to be called.");
    if (idx >= _locs.size()) return false;
    Molecule mol = _mol.lock();
    if (!mol) return false;
    SetAssignment(_results[idx]);
    
    for (auto vs = _g->GetVertices(); vs.first != vs.second; ++vs.first) {
      graph::AGVertex v = *vs.first;
      if (v->IsEdgeMapped()) {
        Bond b = v->GetSourceEdge()->GetBond();
        if (!b) return false;
        switch (v->GetTotalAssigned()) {
          case 2:
            b->SetOrder(BondOrder::SINGLE); break;
          case 4:
            b->SetOrder(BondOrder::DOUBLE); break;
          case 6:
            b->SetOrder(BondOrder::TRIPLE); break;
          case 8:
            b->SetOrder(BondOrder::QUADRUPLE); break;
          default:
            b->SetOrder(BondOrder::UNDEFINED); break;
        }
      } else {
        Atom a = v->GetSourceVertex()->GetAtom();
        if (!a) return false;
        char fc = a->GetElement().GetValenceElectronCount();
        fc -= v->GetTotalAssigned();
        for (auto ns = _g->GetNeighbours(v); ns.first != ns.second; ++ns.first)
          fc -= (*ns.first)->GetTotalAssigned() / 2;
        a->SetFormalCharge(fc);
      }
    }
    return true;
  }
  
  
  void IXElectronAssigner::LoadScoreTable() {
    std::string ddir_path = std::string(IX_DATA_DIRECTORY);
    if (ddir_path.back() != '/') ddir_path.append("/");
    ddir_path.append(settings::ScoreFile);
    utils::FileReader file(ddir_path);
    std::vector<std::string> lines;
    try {
      file.GetAllLines(lines);
    } catch (const std::invalid_argument& e) {
      file.SetFilePath(settings::ScoreFile);
      file.GetAllLines(lines);
    }
    
    const PeriodicTable& PT = GetPeriodicTable();
    
    double multiplier = 100000;
    
    long global_min = std::numeric_limits<long>::max();
    std::map<Element, long> local_atom_min;
    std::map<std::pair<Element, Element>, long> local_bond_min;
    std::set<stdx::quad<Element, Element, char, long>> tmp_dat;
    for (std::string line : lines) {
      std::vector<std::string> items;
      boost::split(items, line, boost::is_any_of(" \t"), boost::token_compress_on);
      
      if (items.size() == 3) {    // Atom score
        // Expect Element, formal charge, score
        const Element& e = PT[items[0]];
        if (!e) continue;
        char fc = std::stoi(items[1]);
        long s = std::llround(std::stod(items[2]) * multiplier);
        // Populate tmp table
        auto tt = stdx::make_quad(e, PT.GetUndefined(), fc, s);
        if (tmp_dat.find(tt) == tmp_dat.end()) tmp_dat.emplace(tt);
        // Determine local/global minimums
        if (local_atom_min.find(e) == local_atom_min.end())
          local_atom_min.emplace(e, s);
        else if (local_atom_min.at(e) > s)
          local_atom_min.at(e) = s;
        if (s < global_min) global_min = s;
        
      } else if (items.size() == 4) {  // Bond score
        // Expect Element A, Element B, Order, score
        Element a = PT[items[0]];
        Element b = PT[items[1]];
        char o = std::stoi(items[2]);
        if (!a || !b || o <= 0) continue;
        long s = std::llround(std::stod(items[3]) * multiplier);
        // Populate tmp table
        auto t1 = stdx::make_quad(a, b, o, s);
        auto t2 = stdx::make_quad(b, a, o, s);
        if (tmp_dat.find(t1) == tmp_dat.end()) tmp_dat.emplace(t1);
        if (tmp_dat.find(t2) == tmp_dat.end()) tmp_dat.emplace(t2);
        // Determine local/global minimums
        std::pair<Element, Element> pa = std::make_pair(a, b);
        std::pair<Element, Element> pb = std::make_pair(b, a);
        if (local_bond_min.find(pa) == local_bond_min.end())
          local_bond_min.emplace(pa, s);
        else if (local_bond_min.at(pa) > s)
          local_bond_min.at(pa) = s;
        if (local_bond_min.find(pb) == local_bond_min.end())
          local_bond_min.emplace(pb, s);
        else if (local_bond_min.at(pb) > s)
          local_bond_min.at(pb) = s;
        if (s < global_min) global_min = s;
      }
    }
    for (auto eets_quad : tmp_dat) {
      if (!(eets_quad.second.GetAtomicNumber())) {
        key_t k = eets_quad.first.GetAtomicNumber();
        k += std::abs(eets_quad.third) << 8;
        if (eets_quad.third < 0) k += 1 << 15;
        _table.emplace(k, eets_quad.fourth - local_atom_min.at(eets_quad.first));
      } else {
        key_t k = eets_quad.first.GetAtomicNumber();
        k += eets_quad.second.GetAtomicNumber() << 8;
        k += (eets_quad.third * 2) << 20;
        auto ee = std::make_pair(eets_quad.first, eets_quad.second);
        _table.emplace(k, eets_quad.fourth - local_bond_min.at(ee));
      }
    }
    _current_file = settings::ScoreFile;
  }
  
  size_t IXElectronAssigner::Run() {
    using op = Option;
    // Check molecule is valid
    Molecule mol = _mol.lock();
    if (!mol) throw std::runtime_error("Molecule is invalid");
    
    auto allow_end = settings::AllowedElements.cend();
    for (auto its = mol->GetAtomIters(); its.first != its.second; ++its.first) {
      Element e = (*its.first)->GetElement();
      if (settings::AllowedElements.find(e) == allow_end)
        throw std::runtime_error("Invalid element in molecule");
    }
    
    // Check settings valid
    static const op valid_pairs = op::Yes|op::No|op::Auto|op::Default;
    static const op valid_cc = op::Yes|op::No|op::Default;
    static const op valid_pre = op::Yes|op::No|op::Default;
    static const op valid_lo_level = op::All|op::Some|op::Default;
    static const op valid_lo_cache = op::All|op::Some|op::None|op::Default;
    if ((settings::ElectronPairs & valid_pairs) != settings::ElectronPairs)
      throw std::runtime_error("Invalid pairs option.");
    if ((settings::ChargedCarbon & valid_cc) != settings::ChargedCarbon)
      throw std::runtime_error("Invalid charged carbon option.");
    if ((settings::Preassign & valid_pre) != settings::Preassign)
      throw std::runtime_error("Invalid pre assign option.");
    if (settings::MaxBondOrder == 0)
      throw std::runtime_error("Maximum bond order of 0 is invalid.");
    
    if (Settings::Algorithm == AssignerAlgorithm::LocalOptimisation) {
      using loset = IXLocalOptimisation::Settings;
      if ((loset::OPTIMISE_LEVEL & valid_lo_level) != loset::OPTIMISE_LEVEL)
        throw std::runtime_error("Invalid optimise level.");
      if ((loset::USE_CACHE & valid_lo_cache) != loset::USE_CACHE)
        throw std::runtime_error("Invalid cache usage.");
    }
    // Load the score table
    if (_table.empty() || settings::ScoreFile != _current_file)
      LoadScoreTable();
    
    // Create the algorithm
    switch (Settings::Algorithm) {
      case AssignerAlgorithm::LocalOptimisation:
        _algo = std::make_unique<IXLocalOptimisation>(_table);
        break;
      case AssignerAlgorithm::AStar:
        _algo = std::make_unique<IXAStarOptimisation>(_table);
        break;
      default:
        throw std::runtime_error("Unsupported optimisation algorithm");
    }
    _algo->Initalise(mol);
    // Run the algorithm
    _algo->Run();
    return _algo->GetOptimalCount();
  }
}


