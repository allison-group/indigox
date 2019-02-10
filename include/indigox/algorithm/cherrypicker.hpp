#include "../classes/athenaeum.hpp"
#include "../classes/forcefield.hpp"
#include "../utils/enum_class_bitwise.hpp"
#include "../utils/fwd_declares.hpp"

#include <cstdint>
#include <list>
#include <type_traits>

#ifndef INDIGOX_ALGORITHM_CHERRYPICKER_HPP
#define INDIGOX_ALGORITHM_CHERRYPICKER_HPP

namespace indigox::algorithm {

  /*! \brief CherryPicker parameterisation algorithm class.
   *  \details The CherryPicker algorithm is a parameterisation algorithm for
   * molecular dynamics simulations. It parameterises novel molecules by
   * identifying fragments of previously parameterised molecules which match
   * portions of the novel molecule, then applying the parameters found to the
   * novel molecule.
   *
   *  Use of CherryPicker is straight forward using the Python bindings, which
   * is the recommended way of doing so. All following code is given as if
   * utilising the Python bindings, and assuming that the following line is
   * included in the script:
   *
   *  \code{.py}
   *  import indigox as ix
   *  \endcode
   *
   * There are two main parts to the CherryPicker algorithm: algorithm setup,
   * and running. Setup consists of initailising a CherryPicker instance,
   * providing an Athenaeum to search through, and changing any of the default
   * settings as desired. The first step of setup requires obtaining a
   * Forcefield which contains all the parameters that could be used to
   * parameterise your target molecule. That is, you are going to parameterise
   * your molecule with this Forcefield. As we primarially use the GROMOS family
   * of forcefields within our group, in particular the <a
   * href="http://dx.doi.org/10.1007/s00249-011-0700-9">54A7 version</a> , a
   * hard coded implementation of the forcefield is provided for simple access.
   * If you wish to use a different version or family of forcefield, consult the
   * Forcefield class for help. With a Forcefield, obtaining a CherryPicker
   * instance then requires simply calling the CherryPicker constructor:
   *
   *  \code{.py}
   *  forcefield = ix.GenerateGROMOS54A7()
   *  cherrypicker = ix.algorithm.CherryPicker(forcefield)
   *  \endcode
   *
   * Left as-is, this CherryPicker instance is useless. It will not be able to
   * parameterise any molecules as it does not have an Athenaeum to search
   * through for matching fragments. As such, the next step is to provide an
   * Athenaeum to CherryPicker. See the Athenaeum class for details on how to
   * create one.
   *
   * \code{.py}
   *  added_ok = cherrypicker.AddAthenaeum(athenaeum)
   *  if not added_ok:
   *    print("Athenaeum failed to add to CherryPicker instance.")
   * \endcode
   *
   *  The AddAthenaeum() method returns a boolean for if the Athenaeum was successfully added to the CherryPicker instance or not. To be successfully added, the Athenaeum::GetForcefield() method must return a Forcefield which matches the Forcefield used to create the CherryPicker instance. As shown above, this return value should be checked to ensure that everything is as expected.
   *
   *  Of further note is that CherryPicker is designed to be able to parameterise molecules using multiple Athenaeum instances.
   */
  class CherryPicker {
  public:
    enum class VertexParameters {
      None = 0,
      ElementType = (1 << 0),
      FormalCharge = (1 << 1),
      CondensedVertices = (1 << 2),
      CyclicNature = (1 << 3),
      Stereochemistry = (1 << 4),
      Aromaticity = (1 << 5),
      Degree = (1 << 6)
    };

    enum class EdgeParameters {
      None = 0,
      BondOrder = (1 << 0),
      Stereochemistry = (1 << 1),
      CyclicNature = (1 << 2),
      Aromaticity = (1 << 3),
      Degree = (1 << 4)
    };

    struct Settings {
      /*! \brief Allows bonds to be parameterised across the overlap region.
       *  \details If one atom of a bond is in the fragment and the other is in
       *  the overlap, setting this to true allows the bond to be parameterised
       *  by the region between the fragment atoms and the overlap atoms. */
      static bool AllowDanglingBonds;

      /*! \brief Allows angles to be parameterised across the overlap region.
       *  \details If two atoms of an angle are in the fragment and the other is
       *  in the overlap, setting this to true allows the angle to be
       *  parameterised by the region between the fragment atoms and the overlap
       *  atoms. */
      static bool AllowDanglingAngles;

      /*! \brief Allows dihedrals to be parameterised across the overlap region.
       *  \details If two atoms of a dihedral are in the fragment and the other
       *  two are in the overlap, setting this to true allows the dihedral to be
       *  parameterised by the region between the fragment atoms and the overlap
       *  atoms. This only applies to dihedrals where the two atoms that are in
       *  the fragment and not the overlap are adjacent to one another. */
      static bool AllowDanglingDihedrals;
      static bool ParameteriseFromAllPermutations;

      static VertexParameters VertexMapping;
      static EdgeParameters EdgeMapping;

      static uint32_t MinimumFragmentSize;
      static uint32_t MaximumFragmentSize;
    };

    CherryPicker() = delete;
    ~CherryPicker() = default;

    CherryPicker(const Forcefield &ff);

    bool AddAthenaeum(const Athenaeum &library);

    bool RemoveAthenaeum(const Athenaeum &library);

    size_t NumAthenaeums() const { return _libs.size(); }

    ParamMolecule ParameteriseMolecule(const Molecule &mol);

    Forcefield GetForcefield() const { return _ff; }

  private:
    //! \brief Forcefield
    Forcefield _ff;
    //! \brief The Athenaeums to use to match
    std::list<Athenaeum> _libs;
  };

  BITWISE_OPERATORS(CherryPicker::VertexParameters);
  BITWISE_OPERATORS(CherryPicker::EdgeParameters);

} // namespace indigox::algorithm

#endif /* INDIGOX_ALGORITHM_CHERRYPICKER_HPP */
