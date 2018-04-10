//
//  electron_optimisation_algorithm.hpp
//  indigox
//
//  Created by Welsh, Ivan on 13/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef ELECTRON_OPTIMISATION_ALGORITHM_HPP
#define ELECTRON_OPTIMISATION_ALGORITHM_HPP

#include <cstdint>
#include <map>
#include <vector>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "../../classes/electron_graph.hpp"
#include "../electron_optimisation.hpp"

namespace indigox {
  namespace algorithm {
    
    typedef boost::dynamic_bitset<> ElnDist;
    
    class ElectronOptimisationAlgorithm {
    protected:
      ElectronOpt* parent_;
      // Result storage
      std::vector<std::pair<ElnDist, uint32_t>> minDistributions_;
      std::map<MolVertPair, ElnVertex> mvp2ev_;
      FCSCORE minScore_;
      FCSCORE upperLimit_;
      
      uint32_t eneBitmask_;   // Used for checking bond energies w/wo charge
      ElnDist previousDist_;
      
    private:
      ElectronOptimisationAlgorithm() = default;
      
    public:
      ElectronOptimisationAlgorithm(ElectronOpt* parent);
      virtual ~ElectronOptimisationAlgorithm() {}
      
      //virtual void Initalise();
      virtual void Run() = 0;
      bool ApplyElectronAssignment(size_t);
      inline size_t GetResultCount() { return minDistributions_.size(); }
      inline FCSCORE GetResultEnergy() { return minScore_; }
      void PopulateMVP2EV();
      
    protected:
      void SetElectronDistribution(ElnDist& dist);
      void DetermineFormalCharges();
      ElnDist CalculateUpperLimit();
      FCSCORE CalculateVertexEnergy(ElnVertex& v);
      FCSCORE CalculateDistributionEnergy(ElnDist);
      
      
    };
    
  }
}

#endif /* ELECTRON_OPTIMISATION_ALGORITHM_HPP */
