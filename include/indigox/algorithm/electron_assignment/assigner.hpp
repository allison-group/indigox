#ifndef INDIGOX_ALGORITHM_ELECTRON_OPTIMISER_HPP
#define INDIGOX_ALGORITHM_ELECTRON_OPTIMISER_HPP

#include <bitset>
#include <map>
#include <set>
#include <vector>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "../../classes/periodictable.hpp"
#include "../../graph/assignment.hpp"
#include "../../utils/common.hpp"
#include "../../utils/numerics.hpp"

// Forward declarations
namespace indigox {
  class IXMolecule;
  using Molecule = std::shared_ptr<IXMolecule>;
  using _Molecule = std::weak_ptr<IXMolecule>;
  namespace algorithm {
    class IXElectronAssigner;
    using ElectronAssigner = std::shared_ptr<IXElectronAssigner>;
  }
}

namespace indigox::algorithm {
  
  //! Type of the score of an electron assignment.
  using score_t = ulong_;
  //! Type of the mask used for assignments.
  using AssignMask = boost::dynamic_bitset<>;
  //! Type of the key for a score table
  using key_t = uint_;
  //! Type of a score table
  using ScoreTable = std::map<key_t, score_t>;
  
  class IXElectronAssigner {
    //! \brief Friendship allows creation of electron assigner instances
    friend ElectronAssigner CreateElectronAssigner(const Molecule& m);
    
  public:
    /*! \brief Base class for an electron assignment algorithm.
     *  \details An AssignAlgorithm is an algorithm for assigning electrons to
     *  an AssignmentGraph. Possible locations to assign electrons are stored
     *  in a vector. If a location can have multiple electrons assigned, it
     *  appears multiple times in the locations vector. An assignment is given
     *  by a bitset with the same length as the locations vector. Each position
     *  in the bitset indicates if the corresponding location should have
     *  electrons assigned, where the position is set if the location should
     *  have electrons assigned, and unset otherwise.
     */
    class AssignAlgorithm {
    protected:
      //! \brief Normal enum for algorithm option indexing.
      enum _bool_opts {
        __pairs,      //!< Use pairs of electrons.
        __charged_c,  //!< Allow charge carbon.
        __preplace,   //!< Preassign electrons.
        __initalised, //!< Algorithm initailised.
        __num_opts    //!< Number of options
      };
      
    public:
      AssignAlgorithm() = delete;   // no default constructor
      
      /*! \brief Normal constructor.
       *  \details Sets the algorithm options from how they are set at the time
       *  of construction.
       *  \param t reference to the score table for the algorithm to utilise. */
      AssignAlgorithm(const ScoreTable& t);
      
      //! \brief Virtual destructor to avoid memory leaks.
      virtual ~AssignAlgorithm() { }
      
      /*! \brief Initalise an assignment algorithm.
       *  \details Initalisation occurs in the following order. First an
       *  AssignmentGraph is generated from the molecule. The number of
       *  electrons to assign are then calculated, accounting for any
       *  preassigned electrons. Finally, the possible locations to assign
       *  electrons are determined. Throughout the initalisation process, sanity
       *  checks are performed to ensure that the initialised state is valid.
       *  \param m the molecule to assign electrons to. */
      virtual void Initalise(const Molecule& m);
      
      //! \brief Required method for algorithms to run.
      virtual void Run() = 0;
      
      /*! \brief Obtain the AssignmentGraph in use by the algorithm.
       *  \details If the molecule assocaited with the assignment algorithm is
       *  no longer valid, or the algorithm has not been initalised, the
       *  returned shared_ptr is empty.
       *  \return the AssignmentGraph used by the algorithm. */
      graph::AssignmentGraph GetAssignmentGraph();
      
      /*! \brief Get the optimised score.
       *  \details If the optimisation has not been run, returned score will be
       *  infinity.
       *  \return the score of the optimised electron assignment. */
      inline score_t GetOptimisedScore() const { return _min_score; }
      
      /*! \brief Get the number of optimal assignments found.
       *  \return the number of optimal assignments found. */
      inline size_ GetOptimalCount() const { return _results.size(); }
      
      /*! \brief Apply an optimised assignment.
       *  \details Applies the optimised assignment at position \p idx to the
       *  Molecule. Application involves determining the bond order and formal
       *  charges from the electron assignment. This method requires the
       *  algorithm the have been initalised.
       *  \param idx the index of the assignment to apply.
       *  \return if the application process was successful or not. */
      bool ApplyAssignment(size_ idx);
      
      /*! \brief If the algorithm has been initalised or not.
       *  \return if the initalise method has been called. */
      bool IsInitalised() const { return _opts[__initalised]; }
      
    protected:
      /*! \brief Calculates the upper score limit.
       *  \details The upper score limit can be used by an assignment algorithm
       *  to improve its execution time. The upper limit is calculated by
       *  iteratively generating an assignment. Iteration involves assigning
       *  electrons to the lowest scoring location, and iterating until all
       *  electrons have been assigned. Once the assignment has been generated,
       *  its score is calculated. If the score is finite, it is incremented by
       *  one. This increment allows easier handling of assignments having a
       *  score the same as the upper limit. This method requires the algorithm
       *  to have been initialised.
       *  \return the assignment that generated the upper limit score.
       *  \throws std::runtime_error if the method is called before initalistion.
       */
      AssignMask CalculateUpperLimit();
      
      /*! \brief Calculate the score of a vertex.
       *  \details Determines the key of the vertex and returns the score of
       *  that key as present in the score table. There are three possibilities
       *  for the score returned. The first is the score present in the score
       *  table. The second is infinity. For a vertex mapped \p v, this occurs
       *  when charged carbons are not allowed and the formal charge on a carbon
       *  is non-zero, when the magnitude of the formal charge on the atom is
       *  larger than that allowed, and if the valency of the atom exceeds that
       *  allowed by its element. For an edge mapped \p v, infinty occurs when
       *  either atom of the bond exceeds its allowed valence. Finally, in
       *  either vertex or edge mapped cases, infinity is returned when the key
       *  is not found in the score table. The final possible score is zero,
       *  which occurs when \p v cannot be assigned any electrons and, in the
       *  vertex mapped case, when none of the neighbours of \p v can have
       *  electrons assigned either. This method requires the algorithm to have
       *  been initialised.
       *  \param v the vertex to calculate the score of.
       *  \return the score of the vertex.
       *  \throws std::runtime_error if the method is called before initalistion.
       */
      score_t CalculateVertexScore(const graph::AGVertex& v) const;
      
      /*! \brief Calculate the score of an assignment.
       *  \details Simply sets the assignment to that provided, and sums all
       *  the scores of the vertices in the AssignmentGraph. If any vertex is
       *  given a score if infinity, the summation halts and infinity is
       *  returned for the score of the assignment. This method requires the
       *  algorithm to have been initialised.
       *  \param a the assignment to calculate the score of.
       *  \return the score of the assignment.
       *  \throws std::runtime_error if the method is called before initalistion.
       */
      score_t CalculateAssignmentScore(const AssignMask& a);
      
      /*! \brief Apply an assignment to the AssignmentGraph.
       *  \details Sets the assigned electron counts for all vertices in the
       *  AssignmentGraph given the provided mask. Only those assignments which
       *  have changed since the previous call to this method are actually
       *  modified. As such, all modification to the AssignmentGraph assigned
       *  counts should only go through this method. This method requires the
       *  algorithm to have been initialised.
       *  \param a the assignment to set.
       *  \throws std::runtime_error if the method is called before initalistion.
       */
      void SetAssignment(const AssignMask& a);
      
    protected:
      //! \brief Molecule working with
      _Molecule _mol;
      //! \brief Assignment graph working on
      graph::AssignmentGraph _g;
      //! \brief Locations to possibly place electrons.
      std::vector<graph::AGVertex> _locs;
      //! \brief Number of electrons to assign.
      size_ _num_e;
      //! \brief State of boolean options
      std::bitset<__num_opts> _opts;
      //! \brief Infinity value
      const score_t _inf;
      //! \brief Score of optimisied assignments
      score_t _min_score;
      //! \brief Upper limit of optimised score
      score_t _limit;
      //! \brief Optimised assignments
      std::vector<AssignMask> _results;
      //! \brief Reference to the score table
      const ScoreTable& _table;
      //! \brief Maximum charge magnitude
      int_ _max_charge;
      //! \brief Maximum number of results
      uint_ _max_results;
      //! \brief Previous assignment mask
      AssignMask _previous_mask;
    };
    
    //! \brief Enum of the various electron assignment optimisation algorithms
    enum class Algorithm {
      LOCAL_OPTIMISATION, //!< The local optimisation method.
      ASTAR,              //!< An A* path finding method.
      FPT                 //!< A dynamic programming method.
    };
    
    struct Settings {
      /*! \brief Which algorithm to use to assign electrons.
       *  \details The default algorithm is the local optimisation method. */
      static Algorithm ALGORITHM;
      
      /*! \brief Assign electrons in pairs instead of singly.
       *  \details The default option (Auto) assigns pairs of electrons when
       *  there is an even number of electrons, and singly if there are an odd
       *  number of electrons. Valid options are Yes, No, Default and Auto. */
      static utils::Option ELECTRON_PAIRS;
      
      /*! \brief Allow carbons to have non zero formal charges.
       *  \details Default option is to allow. Valid options are Yes, No and
       *  Default. */
      static utils::Option CHARGED_CARBON;
      
      /*! \brief Assign some electrons prior to performing optimisation.
       *  \details Default is to do so. Valid options are Yes, No and Default.
       *  Preassigned electrons are fixed and will not moving during the
       *  optimisation. The assigned electrons are as follows:
       *
       *  - Six electrons on one coordinate F, Cl, and Br.
       *  - Four electrons on one coordinate O, and S.
       *  - Two electrons on one coordinate N.
       *  - Four electrons on two coordinate O, and S.
       *  - No electrons on all other types of atoms.
       *  - Two electrons on all bonds. */
      static utils::Option PREASSIGN;
      
      /*! \brief Path to the assignment score file.
       *  \details The assignment score file contains the scores to be assigned
       *  to the various formal charge and bond order states obtained through
       *  an electron assignment. The path should be either relative to the data
       *  directory, or an absolute path. */
      static string_ ASSIGNMENT_SCORE_FILE;
      
      /*! \brief Value of an infinite score.
       *  \details Default value (which probably should not need to be changed)
       *  is std::numeric_limits<score_t>::max(). Any electron assignment
       *  with a final score of the infinite value is considered invalid. */
      static score_t INFINITY_VALUE;
      
      /*! \brief Maximum bond order to assign.
       *  \details This limits the number of electrons which can be assigned
       *  to a bond to twice this value. Default value is 3.
       *  \todo Convert to using the BondOrder enum. */
      static uint_ MAXIMUM_BOND_ORDER;
      
      /*! \brief The maximum allowed charge magnitude on an atom.
       *  \details Any assignment which results in the formal charge of an atom
       *  being larger than this value will be given an infinte score. If this
       *  value is negative, there is no limit applied. Default value is -1. */
      static int_ MAXIMUM_CHARGE_MAGNITUDE;
      
      /*! \brief Maximum number of degenerate score results to calculate.
       *  \details Optimisation algorithms are capable to returning multiple
       *  electron assignments with the same minimum score. This setting
       *  controls the maximum number of such multiples which will be returned.
       *  Small maximum values may result in the optimisation method not giving
       *  a correct minimum. Large values may result in large amounts of memory
       *  and computational time being required. If set to 0, all results will
       *  be returned. Default value is 64. */
      static uint_ MAXIMUM_RESULT_COUNT;
      
      /*! \brief Set of elements for which scores are available.
       *  \details If a molecule contains elements not in this list, the
       *  assignment algorithms will not be able to execute and so an exception
       *  will be thrown when Initalise is called. Only add elements to this
       *  set if they have scores available. The default set of elements, given
       *  the default assignment score file, is: H, C, N, O, S, P, F, Cl, and
       *  Br. */
      static std::set<Element> ALLOWED_ELEMENTS;
    };
    
  public:
    IXElectronAssigner() = delete;  // no default constructor
   
  private:
    /*! \brief Normal constructor. */
    IXElectronAssigner(const Molecule& mol);
    
  public:
    /*! \brief Run the electron assignment.
     *  \details Initalises and runs the assignment algorithm based on the
     *  current settings of both the IXElectronAssigner and the algorithm to run.
     *  Before running, performs a sanity check to ensure that the molecule
     *  meets the allowed elements rules.
     *  \return the number of electron assignments found. */
    size_ Run();
    
    /*! \brief Get the AssignmentGraph used by the assignment algorithm.
     *  \return the AssignmentGraph used by the algorithm. */
    inline graph::AssignmentGraph GetAssignmentGraph() {
      return _algo->GetAssignmentGraph();
    }
    
    /*! \brief Load the current score file.
     *  \details Loads the current set score file. First attempts to load a path
     *  relative to the data directory. If that fails, attempts to load an
     *  absolute path. */
    void LoadScoreTable();
    
    /*! \brief Apply the given assignment.
     *  \details A number of situations will result in unsuscessful application
     *  of an assignment. If the requested index is outside the range of the
     *  available assignments, if the associated molecule is no longer valid, if
     *  any of the expected atom or bond instances are no longer valid. In the
     *  later cases, the assignment may be partially applied.
     *  \param idx the index number of the assignment to apply.
     *  \return if the assignment application was successful or not. */
    inline bool ApplyAssignment(size_ idx) {
      return _algo->ApplyAssignment(idx);
    }
    
  private:
    //! \brief Reference to the currently assigned molecule.
    _Molecule _mol;
    //! \brief Assignment algorithm
    std::unique_ptr<AssignAlgorithm> _algo;
    //! \brief Scores
    ScoreTable _table;
    //! \brief Currently loaded score file
    string_ _current_file;
    
  };
  
  //! \brief Type for the enum of available assignment algorithms
  using AssignerAlgorithm = IXElectronAssigner::Algorithm;
  
  /*! \brief Create an ElectronAssigner.
   *  \param m the molecule this electron assigner is for.
   *  \return a new ElectronAssigner. */
  inline ElectronAssigner CreateElectronAssigner(const Molecule& m) {
    return ElectronAssigner(new IXElectronAssigner(m));
  }
}

#endif /* INDIGOX_ALGORITHM_ELECTRON_OPTIMISER_HPP */
