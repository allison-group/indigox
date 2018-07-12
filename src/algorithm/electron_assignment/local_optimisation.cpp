#include <chrono>

#include <indigox/algorithm/electron_assignment/assigner.hpp>
#include <indigox/algorithm/electron_assignment/local_optimisation.hpp>
#include <indigox/graph/assignment.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox::algorithm {
  
  using Option = utils::Option;
  using AA_Settings = IXElectronAssigner::Settings;
  using LO_Settings = IXLocalOptimisation::Settings;
  
  Option LO_Settings::OPTIMISE_LEVEL = Option::All;
  uint_ LO_Settings::TIMEOUT = 1000;
  Option LO_Settings::USE_CACHE = Option::Default;
  
  IXLocalOptimisation::IXLocalOptimisation(const ScoreTable& t)
  : IXElectronAssigner::AssignAlgorithm(t), _timeout(Settings::TIMEOUT),
  _lo_opts(0) {
    if (Settings::OPTIMISE_LEVEL == Option::All) _lo_opts.set(__all_mins);
    if (Settings::USE_CACHE == Option::All) {
      _lo_opts.set(__cache);
      _lo_opts.set(__cache_inf);
    }
    else if (Settings::USE_CACHE == Option::Some
             || Settings::USE_CACHE == Option::Default) _lo_opts.set(__cache);
  }
  
  void IXLocalOptimisation::Run() {
    using namespace std::chrono;
    high_resolution_clock::time_point start, finish;
    std::map<AssignMask, score_t> scored_assigns;
    std::set<AssignMask> current_min_assigns, round_min_assigns;
    
    AssignMask initial = CalculateUpperLimit();
    score_t current_min_score = _limit;
    if (_limit < _inf) --current_min_score;
    score_t round_min_score = current_min_score;
    current_min_assigns.insert(initial);
    scored_assigns.emplace(initial, current_min_score);
    size_ current_assigns_count = current_min_assigns.size();
    start = high_resolution_clock::now();
    std::vector<AssignMask> nbrs;
    do {
      std::set<AssignMask>::const_iterator begin, end;
      begin = current_min_assigns.begin();
      if (_lo_opts[__all_mins]) end = current_min_assigns.end();
      else end = ++current_min_assigns.begin();
      
      current_min_score = round_min_score;
      round_min_assigns.clear();
      current_assigns_count = current_min_assigns.size();
      
      while (begin != end) {
        if (_timeout > 0) {
          finish = high_resolution_clock::now();
          auto time = duration_cast<milliseconds>(finish - start).count();
          if (time > _timeout) break;
        }
        GetNeighbourAssignments(*begin, nbrs);
        for (AssignMask m : nbrs) {
          score_t score;
          if (_lo_opts[__cache]) {
            auto seen = scored_assigns.find(m);
            if (seen != scored_assigns.end()) score = seen->second;
            else {
              score = CalculateAssignmentScore(m);
              if (score < _inf || _lo_opts[__cache_inf])
                scored_assigns.emplace(m, score);
            }
          } else score = CalculateAssignmentScore(m);
          
          if (score < round_min_score) {
            round_min_assigns.clear();
            round_min_assigns.insert(m);
            round_min_score = score;
          } else if (score == round_min_score) round_min_assigns.insert(m);
        }
        ++begin;
      }
      
      if (round_min_score < current_min_score)
        current_min_assigns.swap(round_min_assigns);
      else if (round_min_score == current_min_score)
        current_min_assigns.insert(round_min_assigns.begin(),
                                   round_min_assigns.end());
      
    } while (current_min_score > round_min_score
             || current_min_assigns.size() != current_assigns_count);
    
    _results.clear();
    _results.reserve(current_min_assigns.size());
    if (current_min_score != _inf) {
      for (AssignMask m : current_min_assigns) {
        _results.emplace_back(m);
        if (_max_results > 0 && _results.size() >= _max_results) break;
      }
    }
    _min_score = current_min_score;
  }
  
  void IXLocalOptimisation::BuildLocationMasks() {
    _loc_masks.clear();
    AssignMask bitmask = AssignMask(_locs.size());
    graph::AGVertex previous;
    for (size_ i = 0; i < _locs.size(); ++i) {
      graph::AGVertex current = _locs[i];
      if (current == previous) bitmask.set(i);
      else {
        if (i) _loc_masks.emplace(previous, bitmask);
        bitmask.reset();
        bitmask.set(i);
        previous = current;
      }
    }
    if (bitmask.count()) _loc_masks.emplace(previous, bitmask);
  }
  
  void IXLocalOptimisation::GetNeighbourAssignments(const AssignMask &m,
                                                std::vector<AssignMask> &nbrs) {
    nbrs.clear();
    LocMasks_t::const_iterator srcBegin = _loc_masks.begin();
    LocMasks_t::const_iterator srcEnd = _loc_masks.end();
    
    for (; srcBegin != srcEnd; ++srcBegin) {
      size_ srcIdx = srcBegin->second.find_first();
      size_ srcLoc = srcIdx + (m & srcBegin->second).count() - 1;
      if (srcIdx > srcLoc) continue;
      LocMasks_t::const_iterator tgtBegin = _loc_masks.begin();
      LocMasks_t::const_iterator tgtEnd = _loc_masks.end();
      for (; tgtBegin != tgtEnd; ++tgtBegin) {
        if (tgtBegin == srcBegin) continue;
        size_ tgtIdx = tgtBegin->second.find_first();
        size_ tgtCount = tgtBegin->second.count();
        size_ tgtLoc = tgtIdx + (m & tgtBegin->second).count();
        if (tgtCount + tgtIdx == tgtLoc) continue;
        AssignMask nbr(m);
        nbr.flip(srcLoc);
        nbr.flip(tgtLoc);
        nbrs.emplace_back(nbr);
      }
    }
  }
  
  void IXLocalOptimisation::Initalise(const indigox::Molecule &mol) {
    IXElectronAssigner::AssignAlgorithm::Initalise(mol);
    BuildLocationMasks();
  }
  
}
