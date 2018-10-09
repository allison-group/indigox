#include <cstdint>
#include <list>

#include "../utils/fwd_declares.hpp"

#ifndef INDIGOX_ALGORITHM_CHERRYPICKER_HPP
#define INDIGOX_ALGORITHM_CHERRYPICKER_HPP

namespace indigox::algorithm {
  
  class IXCherryPicker {
  public:
    
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
    };
    
    IXCherryPicker() = delete;
    ~IXCherryPicker() = default;
    
    IXCherryPicker(Forcefield ff) : _ff(ff) { }
    
    bool AddAthenaeum(Athenaeum library);
    
    bool RemoveAthenaeum(Athenaeum library);
    
    size_t NumAthenaeums() const { return _libs.size(); }
    
    ParamMolecule ParameteriseMolecule(Molecule mol) const;
    
    Forcefield GetForcefield() const { return _ff; }
  private:
    //! \brief Forcefield
    Forcefield _ff;
    //! \brief The Athenaeums to use to match
    std::list<Athenaeum> _libs;
  };
  
}

#endif /* INDIGOX_ALGORITHM_CHERRYPICKER_HPP */
