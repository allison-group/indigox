//
//  electron_optimisation.hpp
//  indigox
//
//  Created by Welsh, Ivan on 12/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef ELECTRON_OPTIMISATION_HPP
#define ELECTRON_OPTIMISATION_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include "../classes/electron_graph.hpp"
#include "../classes/molecular_graph.hpp"

namespace indigox {
  namespace algorithm {
    class ElectronOptimisationAlgorithm;
    class LocalOptimisation;
    class AStarOptimisation;
    class FPTOptimisation;
  }
  
  class ElectronOpt {
    friend class indigox::algorithm::ElectronOptimisationAlgorithm;
    friend class indigox::algorithm::LocalOptimisation;
    friend class indigox::algorithm::AStarOptimisation;
    friend class indigox::algorithm::FPTOptimisation;
    
  public:
    uint32_t electronsToAdd_;
    std::vector<MolVertPair> possibleLocations_;
    MolecularGraph molGraph_;
    ElectronGraph elnGraph_;
    std::unordered_map<uint32_t, FCSCORE> scores_;
    std::shared_ptr<algorithm::ElectronOptimisationAlgorithm> algo_;
    FCSCORE finalScore_;
    
  public:
    ElectronOpt();
    ElectronOpt(MolecularGraph G);
    
  public:
    void SetMolecularGraph(MolecularGraph G);
    size_t Run();
    inline FCSCORE GetMinimisedEnergy() { return finalScore_; }
    bool ApplyElectronAssigment(size_t);
    
  private:
    void DetermineElectronsToAdd();
    void DeterminePotentialElectronLocations();
    void LoadScores();
    void SortPotentialLocations();
    
  };
}

#endif /* ELECTRON_OPTIMISATION_HPP */
