/** @file options.hpp
 *  @brief Declaration of all available options
 *  @author Ivan Welsh
 *  @date 7 January 2018
 *  @lastmodify 7 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#include <cstdint>
#include <string>
#include <set>

#ifndef INDIGOX_UTILS_OPTIONS_HPP
#define INDIGOX_UTILS_OPTIONS_HPP

#ifndef FCSCORE
#define FCSCORE uint32_t
#endif

namespace indigox {
  
  /** @struct Options options.hpp utils/options.hpp
   *  @brief Static struct containing all runtime changeable options for
   *  indigoX.
   *  @since 0.1
   */
  struct Options {
    
    struct AssignElectrons {
      struct AStar {
        enum class Heuristic {
          PROMISCUOUS,
          ABSTEMIOUS
        };
        static Heuristic HEURISTIC;
        static uint64_t MEGABYTE_LIMIT;
      };
      
      struct FPT {
        enum class PermAlgo {
          RANDOM,
          QUICKBB,
          MINDEGREE,
          MINADDEDGES
        };
        static std::string LIBTW_JAR_FILE;
        static bool ADD_EDGES_TO_TD;
        static PermAlgo PERM_ALGO;
        static uint32_t MINIMUM_PROPAGATION_DEPTH;
      };
      
      struct LocalOptimisation {
        static bool OPTIMISE_ALL_MINIMUMS;
        static bool CACHE_RESULTS;
        static bool CACHE_INFINITIES;
        static uint32_t TIMEOUT_LIMIT;  // milliseconds
      };
      
      enum class Algorithm {
        LOCAL_OPTIMISATION,
        ASTAR,
        FPT
      };
      
      static Algorithm ALGORITHM;
      
      static std::string ATOM_ENERGY_FILE;
      static std::string BOND_ENERGY_FILE;
      
      static FCSCORE INF;
      static uint8_t MAXIMUM_BOND_ORDER;
      static bool USE_ELECTRON_PAIRS;
      static bool AUTO_USE_ELECTRON_PAIRS;
      static bool USE_CHARGED_BOND_ENERGIES;
      static uint32_t HIGHEST_MAGNITUDE_CHARGE;
      static bool ALLOW_CHARGED_CARBON;
      static std::set<std::string> ALLOWED_ELEMENTS;
      static uint32_t MAXIMUM_RESULT_COUNT;
      static bool PREPLACE_ELECTRONS;
    };
    
    
    
    /// @brief Default data directory path.
    static std::string DATA_DIRECTORY;
    
    /// @brief Reset all attributes to their default values.
    static void Reset();
  };
  
}

#endif /* INDIGOX_UTILS_OPTIONS_HPP */
