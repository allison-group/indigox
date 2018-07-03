#include <bitset>
#include <map>
#include <vector>

#include "assigner.hpp"
#include "../../graph/assignment.hpp"
#include "../../utils/common.hpp"
#include "../../utils/numerics.hpp"

#ifndef INDIGOX_ALGORITHM_ELECTRON_ASSIGNMENT_LOCAL_OPTIMISATION_HPP
#define INDIGOX_ALGORITHM_ELECTRON_ASSIGNMENT_LOCAL_OPTIMISATION_HPP

namespace indigox::algorithm {
  
  /*! \brief Local optimisation electron assignment algorithm.
   *  \details Stuff
   */
  class IXLocalOptimisation : public IXElectronAssigner::AssignAlgorithm {
    //! \brief Normal enum for local optimisation algorithm option indexing.
    enum _lo_bool_opts {
      __all_mins,   //!< Minimise from all assignments
      __cache,      //!< Cache previously seen assignments
      __cache_inf,  //!< Cache previously seen assignments with infinte score
      __num_opts    //!< Number of options
    };
  public:
    struct Settings {
      /*! \brief Amount of optimisation to perform.
       *  \details At each round to optimisation, how many of the current
       *  minimum score assignments should be optimised down. Valis options are
       *  All, Some and Default. Default is equivalent to Some. With All, all
       *  of the minimum score assignments are optimised from. With Some, only
       *  the first minimum score assignment is optimised down. */
      static utils::Option OPTIMISE_LEVEL;
      
      /*! \brief Time limit for the optimisation, in milliseconds.
       *  \details If the optimisation takes longer than this amount of time, it
       *  is halted and whatever the current minimum score is is set as the
       *  optimised score. If set to 0, no timeout is applied. Default value is
       *  5000. */
      static uint_ TIMEOUT;
      
      /*! \brief Amount of cache to use.
       *  \details Cache the scores of seen assignments so that they do not need
       *  to be recalculated as they reoccur. Valid options are All, Some, None,
       *  and Default. Default is equivalent to Some. With All, all seen
       *  assignments are cached. With Some, only those assignments that are
       *  not infinite scores are cached. With None, no caching occurs. */
      static utils::Option USE_CACHE;
    };
    
  public:
    IXLocalOptimisation() = delete;  // no default constructor
    
    /*! \brief Normal constructor.
     *  \details Sets the options of the algorithm from their state at the time
     *  of construction. Also calls the base constructor.
     *  \param t the score table to use. */
    IXLocalOptimisation(const ScoreTable& t);
    
    /*! \brief Initalisation method.
     *  \details Overrides the base class initalisation method so that location
     *  masks can be generated as part of initalisation.
     *  \param mol the Molecule to initalise with. */
    virtual void Initalise(const Molecule& mol) override;
    
    
    /*! \brief Run the algorithm.
     *  \details Runs the algorithm and populates the results. */
    virtual void Run() override;
  
  private:
    /*! \brief Build the location masks.
     *  \details Location masks are used to indicate the positions of the
     *  possible electron location vertices within the \p _locs vector. */
    void BuildLocationMasks();
    
    /*! \brief Get the neighbouring assignments of an assignment.
     *  \details A neghbour of an assignment is an assignment which can be
     *  obtained by unassigning electrons from one location and assigning them
     *  to another location. The neighbours of the assignment are all such
     *  assignments that can be created from a given assignment.
     *  \param m the assignment to get the neighbours of.
     *  \param[out] nbrs vector to store the neighbours in. */
    void GetNeighbourAssignments(const AssignMask& m,
                                 std::vector<AssignMask>& nbrs);
  
    
  private:
    //! \brief Type to use to store location masks
    using LocMasks_t = std::map<graph::AGVertex, AssignMask>;
    //! \brief Location mask storage.
    LocMasks_t _loc_masks;
    //! \brief Timeout limit
    uint_ _timeout;
    //! \brief LocalOptimsation options
    std::bitset<__num_opts> _lo_opts;
  };
  
}

#endif /* INDIGOX_ALGORITHM_ELECTRON_ASSIGNMENT_LOCAL_OPTIMISATION_HPP */
