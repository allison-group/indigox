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

    size_t NumAthenaeums() const {
      return _libs.size();
    }

    ParamMolecule ParameteriseMolecule(const Molecule &mol);

    Forcefield GetForcefield() const {
      return _ff;
    }

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
