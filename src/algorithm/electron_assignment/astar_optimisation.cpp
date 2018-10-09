#include <algorithm>
#include <array>
#include <cstdint>
#include <queue>
#include <string>
#include <vector>

#include <boost/algorithm/string/join.hpp>

#include <indigox/algorithm/electron_assignment/assigner.hpp>
#include <indigox/algorithm/electron_assignment/astar_optimisation.hpp>
#include <indigox/classes/atom.hpp>
#include <indigox/classes/bond.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/graph/molecular.hpp>

namespace indigox::algorithm {
  
  std::string __PrintVertexName(const graph::AGVertex& v) {
    if (v->IsVertexMapped())
      return v->GetSourceVertex()->GetAtom()->GetName();
    std::vector<std::string> names; names.reserve(2);
    names.emplace_back(v->GetSourceEdge()->GetBond()->GetSourceAtom()->GetName());
    names.emplace_back(v->GetSourceEdge()->GetBond()->GetTargetAtom()->GetName());
    return boost::join(names, ", ");
  }
  
  using AA_settings = IXElectronAssigner::Settings;
  using AS_settings = IXAStarOptimisation::Settings;
  using Heuristic = IXAStarOptimisation::Heuristic;
  
  uint64_t AS_settings::MEMORY_LIMIT = 1024;
  Heuristic AS_settings::HEURISTIC = Heuristic::Abstemious;
  
  enum class SortOrder {
    F_ONE,
    CL_ONE,
    BR_ONE,
    O_ONE,
    S_ONE,
    O_TWO,
    S_TWO,
    C_TWO_N_ONE,
    C_TWO_C_TWO,
    C_THREE_C_THREE,
    N_ONE,
    P_ONE,
    N_TWO,
    P_TWO,
    N_THREE,
    P_THREE,
    C_THREE_O_ONE,
    O_ONE_S_FOUR,
    C_THREE_S_ONE,
    O_ONE_P_FOUR,
    C_THREE_N_TWO,
    N_TWO_N_TWO,
    N_THREE_O_ONE,
    UNDEFINED
  };
  
  SortOrder GetVertexSortOrder(const graph::AGVertex &v) {
    if (v->IsVertexMapped()) {
      Atom a = v->GetSourceVertex()->GetAtom();
      switch (a->GetElement()->GetAtomicNumber()) {
        case 7:   // Nitrogen
          switch (a->NumBonds()) {
            case 1: return SortOrder::N_ONE;
            case 2: return SortOrder::N_TWO;
            case 3: return SortOrder::N_THREE;
            default: return SortOrder::UNDEFINED;
          }
        case 8:   // Oxygen
          switch (a->NumBonds()) {
            case 1: return SortOrder::O_ONE;
            case 2: return SortOrder::O_TWO;
            default: return SortOrder::UNDEFINED;
          }
        case 9:   // Fluorine
          switch (a->NumBonds()) {
            case 1: return SortOrder::F_ONE;
            default: return SortOrder::UNDEFINED;
          }
        case 15:  // Phosphorus
          switch (a->NumBonds()) {
            case 1: return SortOrder::P_ONE;
            case 2: return SortOrder::P_TWO;
            case 3: return SortOrder::P_THREE;
            default: return SortOrder::UNDEFINED;
          }
        case 16:  // Sulfur
          switch (a->NumBonds()) {
            case 1: return SortOrder::S_ONE;
            case 2: return SortOrder::S_ONE;
            default: return SortOrder::UNDEFINED;
          }
        case 17:  // Chlorine
          switch (a->NumBonds()) {
            case 1: return SortOrder::CL_ONE;
            default: return SortOrder::UNDEFINED;
          }
        case 35:  // Bromine
          switch (a->NumBonds()) {
            case 1: return SortOrder::BR_ONE;
            default: return SortOrder::UNDEFINED;
          }
        default: return SortOrder::UNDEFINED;
      }
    } else {
      Atom a = v->GetSourceEdge()->GetBond()->GetSourceAtom();
      Atom b = v->GetSourceEdge()->GetBond()->GetTargetAtom();
      switch (a->GetElement()->GetAtomicNumber()) {
        case 6:
          switch (b->GetElement()->GetAtomicNumber()) {
            case 6:
              if (a->NumBonds() == 2 && b->NumBonds() == 2)
                return SortOrder::C_TWO_C_TWO;
              else if (a->NumBonds() == 3 && b->NumBonds() == 3)
                return SortOrder::C_THREE_C_THREE;
              else return SortOrder::UNDEFINED;
            case 7:
              if (a->NumBonds() == 2 && b->NumBonds() == 1)
                return SortOrder::C_TWO_N_ONE;
              else if (a->NumBonds() == 3 && b->NumBonds() == 2)
                return SortOrder::C_THREE_N_TWO;
              else return SortOrder::UNDEFINED;
            case 8:
              if (a->NumBonds() == 3 && b->NumBonds() == 1)
                return SortOrder::C_THREE_O_ONE;
              else return SortOrder::UNDEFINED;
            case 16:
              if (a->NumBonds() == 3 && b->NumBonds() == 1)
                return SortOrder::C_THREE_S_ONE;
              else return SortOrder::UNDEFINED;
            default: return SortOrder::UNDEFINED;
          }
        case 7:
          switch (b->GetElement()->GetAtomicNumber()) {
            case 6:
              if (a->NumBonds() == 1 && b->NumBonds() == 2)
                return SortOrder::C_TWO_N_ONE;
              else if (a->NumBonds() == 2 && b->NumBonds() == 3)
                return SortOrder::C_THREE_N_TWO;
              else return SortOrder::UNDEFINED;
            case 7:
              if (a->NumBonds() == 2 && b->NumBonds() == 2)
                return SortOrder::N_TWO_N_TWO;
              else return SortOrder::UNDEFINED;
            case 8:
              if (a->NumBonds() == 3 && b->NumBonds() == 1)
                return SortOrder::N_THREE_O_ONE;
              else return SortOrder::UNDEFINED;
            default: return SortOrder::UNDEFINED;
          }
        case 8:
          switch (b->GetElement()->GetAtomicNumber()) {
            case 6:
              if (a->NumBonds() == 1 && b->NumBonds() == 3)
                return SortOrder::C_THREE_O_ONE;
              else return SortOrder::UNDEFINED;
            case 7:
              if (a->NumBonds() == 1 && b->NumBonds() == 3)
                return SortOrder::N_THREE_O_ONE;
              else return SortOrder::UNDEFINED;
            case 15:
              if (a->NumBonds() == 1 && b->NumBonds() == 4)
                return SortOrder::O_ONE_P_FOUR;
              else return SortOrder::UNDEFINED;
            case 16:
              if (a->NumBonds() == 1 && b->NumBonds() == 4)
                return SortOrder::O_ONE_S_FOUR;
              else return SortOrder::UNDEFINED;
            default: return SortOrder::UNDEFINED;
          }
        case 15:
          switch (b->GetElement()->GetAtomicNumber()) {
            case 8:
              if (a->NumBonds() == 4 && b->NumBonds() == 1)
                return SortOrder::O_ONE_P_FOUR;
              else return SortOrder::UNDEFINED;
            default: return SortOrder::UNDEFINED;
          }
        case 16:
          switch (b->GetElement()->GetAtomicNumber()) {
            case 6:
              if (a->NumBonds() == 1 && b->NumBonds() == 3)
                return SortOrder::C_THREE_S_ONE;
              else return SortOrder::UNDEFINED;
            case 8:
              if (a->NumBonds() == 4 && b->NumBonds() == 1)
                return SortOrder::O_ONE_S_FOUR;
              else return SortOrder::UNDEFINED;
            default: return SortOrder::UNDEFINED;
          }
        default: return SortOrder::UNDEFINED;
      }
    }
  }

  bool QueueItem::operator>(const QueueItem &r) const  {
    if (IsInfinite() && r.IsInfinite()) return false;
    else if (IsInfinite()) return true;
    else if (r.IsInfinite()) return false;
    else if (Total() == r.Total()) return CalcCount() > r.CalcCount();
    else return Total() > r.Total();
  }
  
  IXAStarOptimisation::IXAStarOptimisation(const ScoreTable& t)
  : IXElectronAssigner::AssignAlgorithm(t) { }
  
  void IXAStarOptimisation::Run() {
    NbrAssigns nbrs;
    CalculateUpperLimit();
    score_t target_score = _limit;
    for (graph::AGVertex v : _locs) std::cout << __PrintVertexName(v) << std::endl;
    while (!_q.empty()) {
      QueueItem source = _q.top();
      _q.pop();
      
      // Check exit conditions
      if (source.unchange_mask.all()) {
        if (target_score == _limit) target_score = source.path;
        if (source.path == target_score) {
          _results.emplace_back(source.assignment);
          if (_max_results && _results.size() >= _max_results) break;
        } else if (source.Total() > target_score) break;
        continue;
      } else if (target_score != _limit && source.path > target_score) break;
      
      // Generate the next queue items
      std::cout << __PrintVertexName(_locs[source.nbr_begin_idx]) << std::endl;
      size_t num_nbrs = GenerateNextAssignments(source, nbrs);
      for (size_t i = 0; i < num_nbrs; ++i) {
        graph::AGVertex nbr_vertex = _locs[source.nbr_begin_idx];
        size_t nbr_index = source.nbr_begin_idx + _loc_counts[GetLocationPosition(nbr_vertex)];
        QueueItem item(source.path, source.heuristic, nbrs[i],
                       source.unchange_mask, source.calc_mask,
                       LocMask(source.calc_mask.size()), nbr_index);
        // Unchange_mask and calc_mask are going to be the same between neighbours
        // So should be able to move outside of loop for very minor speed up
        if (item.assignment.count() == _num_e) {
          item.unchange_mask.set();
          item.calc_mask.set();
          item.nbr_begin_idx = _locs.size();
        } else {
          item.unchange_mask.set(GetLocationPosition(nbr_vertex));
          item.calc_mask = DetermineCalculableLocations(item);
        }
        item.new_calc_mask = item.calc_mask - source.calc_mask;
        score_t new_path = CalculateNewPathScore(item);
        if (new_path != _inf) item.path += new_path;
        else item.path = _inf;
        item.heuristic = CalculateHeuristicScore(item);
        
        if (item.path < _inf && item.heuristic < _inf
            && (item.path + item.heuristic) < _limit)
          _q.emplace(item);
      }
      if (_len_limit > 0 && _q.size() > _len_limit) break;
    }
    _q.clear();
    if (target_score != _limit) _min_score = target_score;
    else _min_score = _inf;
  }
  
  void IXAStarOptimisation::DetermineLocationCounts() {
    auto verts = _g->GetVertices();
    _unique_locs.insert(_unique_locs.begin(), verts.first, verts.second);

    //  Need to order the unique locations such that the vertex mapped ones
    //  come first.
    auto loc_order = [](const graph::AGVertex& a, const graph::AGVertex& b) {
      if (a->IsVertexMapped() && b->IsEdgeMapped()) return true;
      else return false;
    };
    std::sort(_unique_locs.begin(), _unique_locs.end(), loc_order);
    
    _loc_counts.reserve(_unique_locs.size());
    for (graph::AGVertex v : _unique_locs)
      _loc_counts.emplace_back(std::count(_locs.begin(), _locs.end(), v));
  }
  
  void IXAStarOptimisation::DetermineRequiredUnchangeables() {
    size_t num_v = _g->NumVertices();
    _req_unchange.reserve(num_v);
    
    for (size_t pos = 0; pos < num_v; ++pos) {
      graph::AGVertex v = _unique_locs[pos];
      _req_unchange.emplace_back(num_v);
      auto nbrs = _g->GetNeighbours(v);
      for (; nbrs.first != nbrs.second; ++nbrs.first)
        _req_unchange.back().set(GetLocationPosition(*nbrs.first));
    }
  }
  
  void IXAStarOptimisation::DetermineInitialAssignment() {
    QueueItem item(_locs.size(), _unique_locs.size());
    item.unchange_mask.set();
    for (size_t i = 0; i < _loc_counts.size(); ++i) {
      if (_loc_counts[i]) item.unchange_mask.reset(i);
    }
    
    item.calc_mask = DetermineCalculableLocations(item);
    item.new_calc_mask = item.calc_mask;
    
    item.path = CalculateNewPathScore(item);
    item.heuristic = CalculateHeuristicScore(item);
    _q.emplace(item);
  }
  
  LocMask IXAStarOptimisation::DetermineCalculableLocations(const QueueItem &q) const {
    LocMask calculable(q.calc_mask);
    if (q.unchange_mask.all() || q.assignment.count() >= _num_e) {
      calculable.set();
      return calculable;
    }
    
    size_t pos = q.unchange_mask.find_first();
    while (pos < q.unchange_mask.size()) {
      if (!calculable.test(pos)) {
        if ((q.unchange_mask & _req_unchange[pos]) == _req_unchange[pos])
          calculable.set(pos);
      }
      pos = q.unchange_mask.find_next(pos);
    }
    return calculable;
  }
  
  score_t IXAStarOptimisation::CalculateNewPathScore(const QueueItem &q) {
    if (_num_e < q.assignment.count()) return _inf;
    if ((q.assignment.size() - q.nbr_begin_idx)
        < (_num_e - q.assignment.count())) return _inf;
      
    SetAssignment(q.assignment);
    score_t tot_score = 0;
    size_t pos = q.new_calc_mask.find_first();
    while (pos < q.new_calc_mask.size()) {
      score_t score = CalculateVertexScore(_unique_locs[pos]);
      if (score == _inf) return _inf;
      tot_score += score;
      pos = q.new_calc_mask.find_next(pos);
    }
    return tot_score;
  }
  
  score_t IXAStarOptimisation::CalculateHeuristicScore(const QueueItem &q) const {
    if (q.assignment.count() == _num_e) return 0;
    if (q.assignment.count() > _num_e) return _inf;
    if (Settings::HEURISTIC == Heuristic::Promiscuous)
      return Promiscuous(q);
    if (Settings::HEURISTIC == Heuristic::Abstemious)
      return Abstemious(q);
    return 0;
  }
  
  score_t IXAStarOptimisation::Promiscuous(const QueueItem &q) const {
    score_t h = 0;
    LocMask uncalcuable(q.calc_mask);
    uncalcuable.flip();
    size_t pos = uncalcuable.find_first();
    while (pos < uncalcuable.size()) {
      graph::AGVertex v = _unique_locs[pos];
      score_t min_s = _inf;
      if (v->IsVertexMapped()) {
        key_t k1 = v->GetSourceVertex()->GetAtom()->GetElement()->GetAtomicNumber();
        for (key_t fc = 0; fc < 10; ++fc) {
          // Positive formal charges
          key_t mask = k1 + (fc << 8);
          if (_table.find(mask) != _table.end()) {
            score_t s = _table.at(mask);
            if (s < min_s) min_s = s;
          }
          // Negative formal charges
          mask += (1 << 15);
          if (_table.find(mask) != _table.end()) {
            score_t s = _table.at(mask);
            if (s < min_s) min_s = s;
          }
        }
      } else {
        auto nbrs = _g->GetNeighbours(v);
        graph::AGVertex u1 = *nbrs.first++;
        graph::AGVertex u2 = *nbrs.first;
        key_t k1 = u1->GetSourceVertex()->GetAtom()->GetElement()->GetAtomicNumber();
        k1 += (u2->GetSourceVertex()->GetAtom()->GetElement()->GetAtomicNumber() << 8);
        for (key_t bond_e = 1; bond_e < 10; ++bond_e) {
          key_t mask = k1 + (bond_e << 20);
          if (_table.find(mask) == _table.end()) continue;
          score_t s = _table.at(mask);
          if (s < min_s) min_s = s;
        }
      }
      if (min_s != _inf) h += min_s;
      else return _inf;
      pos = uncalcuable.find_next(pos);
    }
    return h;
  }
  
  score_t IXAStarOptimisation::Abstemious(const QueueItem &q) const {
    // Note where electrons can be placed
    size_t extras = _num_e - q.assignment.count();
    std::vector<uint8_t> extra_positions;
    extra_positions.reserve(_unique_locs.size());
    for (size_t i = 0; i < _unique_locs.size(); ++i) {
      if (q.unchange_mask[i]) extra_positions.emplace_back(0);
      else if (extras <= _loc_counts[i]) extra_positions.emplace_back(extras);
      else extra_positions.emplace_back(_loc_counts[i]);
    }
    
    if (_opts[__pairs]) extras += extras;
    score_t h = 0;
    LocMask uncalcuable(q.calc_mask);
    uncalcuable.flip();
    
    size_t pos = uncalcuable.find_first();
    while (pos < uncalcuable.size()) {
      graph::AGVertex v = _unique_locs[pos];
      score_t min_s = _inf;
      if (v->IsVertexMapped()) {
        size_t nbr_e = 0;
        size_t me_e = extra_positions[pos];
        for (auto n = _g->GetNeighbours(v); n.first != n.second; ++n.first)
          nbr_e += extra_positions[GetLocationPosition(*n.first)];
        if (_opts[__pairs]) me_e += me_e;
        else nbr_e /= 2;
        size_t step = 1;
        if (!(me_e % 2) && ! nbr_e) step = 2;
        size_t toplace = nbr_e + me_e;
        if (extras < toplace) toplace = extras;
        Element e = v->GetSourceVertex()->GetAtom()->GetElement();
        int fc = e->GetValenceElectronCount() - v->GetPreAssignedCount();
        key_t k1 = e->GetAtomicNumber();
        for (int i = 0; i <= static_cast<int>(toplace); i += step) {
          key_t mask = k1 + (abs(fc + i) << 8);
          if ((fc + i) < 0) mask += (1 << 15);
          if (_table.find(mask) == _table.end()) continue;
          score_t s = _table.at(mask);
          if (s < min_s) min_s = s;
        }
        
      } else {
        auto nbrs = _g->GetNeighbours(v);
        graph::AGVertex u1 = *nbrs.first++;
        graph::AGVertex u2 = *nbrs.first;
        key_t k1 = u1->GetSourceVertex()->GetAtom()->GetElement()->GetAtomicNumber();
        k1 += (u2->GetSourceVertex()->GetAtom()->GetElement()->GetAtomicNumber() << 8);
        size_t toplace = extra_positions[pos];
        if (_opts[__pairs]) toplace += toplace;
        
        for (size_t e = 0; e <= toplace; ++e) {
          key_t mask = k1 + ((v->GetPreAssignedCount() + e) << 20);
          if (_table.find(mask) == _table.end()) continue;
          score_t s = _table.at(mask);
          if (s < min_s) min_s = s;
        }
      }
      if (min_s < _inf) h += min_s;
      else return _inf;
      pos = uncalcuable.find_next(pos);
    }
    return h;
  }
  
  size_t IXAStarOptimisation::GenerateNextAssignments(const QueueItem &q,
                                                     NbrAssigns &out) const {
    graph::AGVertex nbr_v = _locs[q.nbr_begin_idx];
    size_t nbr_count = _loc_counts[GetLocationPosition(nbr_v)] + 1;
    for (size_t i = 0; i < nbr_count; ++i) {
      AssignMask nbr(q.assignment);
      for (size_t j = 0; j < i; ++j) nbr.set(j + q.nbr_begin_idx);
      out[i] = nbr;
    }
    return nbr_count;
  }
  
  
  void IXAStarOptimisation::Initalise(const Molecule &mol) {
    IXElectronAssigner::AssignAlgorithm::Initalise(mol);
    
    // Sort the possible locations
    auto cmp = [](const graph::AGVertex& l, const graph::AGVertex& r) {
      return GetVertexSortOrder(l) < GetVertexSortOrder(r);
    };
    std::stable_sort(_locs.begin(), _locs.end(), cmp);
    
//    for (graph::AGVertex v : _locs) std::cout << __PrintVertexName(v) << std::endl;
    
    // Setup all the containers and stuff
    DetermineLocationCounts();
    DetermineRequiredUnchangeables();
    std::cout << "Required Unchangables: " << std::endl;
    for (size_t i = 0; i < _unique_locs.size(); ++i) {
      std::cout << __PrintVertexName(_unique_locs[i]) << ": " << _req_unchange[i] << std::endl;
    }
    
    // Figure out the _len_limit
    if (Settings::MEMORY_LIMIT == 0) _len_limit = 0;
    else {
      QueueItem t;
      t.assignment = AssignMask(_locs.size());
      t.calc_mask = LocMask(_g->NumVertices());
      t.unchange_mask = LocMask(_g->NumVertices());
      
      size_t per_item = sizeof(QueueItem);
      per_item += sizeof(AssignMask::block_type) * t.assignment.num_blocks();
      per_item += sizeof(LocMask::block_type) * t.calc_mask.num_blocks();
      per_item += sizeof(LocMask::block_type) * t.unchange_mask.num_blocks();
      
      size_t count = 1024 * 1024 / per_item;
      _len_limit = Settings::MEMORY_LIMIT * count;
      _q.reserve(_len_limit);
    }
    
    // Create initial assignment
    DetermineInitialAssignment();
    
  }
  
  
}
