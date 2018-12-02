/*! \file angle.hpp */
#include <array>
#include <memory>
#include <string>

#include "forcefield.hpp"
#include "../utils/common.hpp"
#include "../utils/fwd_declares.hpp"
#include "../utils/counter.hpp"
#include "../utils/triple.hpp"

#ifndef INDIGOX_CLASSES_ANGLE_HPP
#define INDIGOX_CLASSES_ANGLE_HPP

namespace indigox {
  class Angle {
    //! \brief Friendship allows Molecule to create new angles.
    friend class indigox::Molecule;
    //! \brief Friendship allows serialisation.
    friend class cereal::access;
    
  public:
    //! \brief Container for storing Atom reference assigned to an Angle.
    using AngleAtoms = std::array<Atom, 3>;
    
  private:
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Angle);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Angle, ang);
    
  private:
    /*! \brief Normal constructor.
     *  \details Creates an angle between three atoms, linking it to the given
     *  Molecule. Atom \p b is always the central atom of the angle.
     *  \param a,b,c the three atoms to construct an angle between.
     *  \param m the molecule to assign the angle to. */
    Angle(const Atom& a, const Atom& b, const Atom& c, const Molecule& m);
    
  public:
    /*! \brief Tag of the angle.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the angle, use the Angle::GetUniqueID()
     *  method.
     *  \return the tag assigned to the angle. */
    int32_t GetTag() const;
    
    /*! \brief Molecule this angle is associated with.
     *  \details The returned shared_ptr is empty if the angle is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this angle. */
    const Molecule& GetMolecule() const;
    
    /*! \brief String representation of the angle.
     *  \details If any of the atoms of the angle are not valid, the returned
     *  string is Angle(MALFORMED), otherwise it is of the form:
     *  Angle(NAME_BEGIN, NAME_CENTRAL, NAME_END).
     *  \return a string representation of the dihedral. */
    std::string ToString() const;
    
    /*! \brief Set the tag of this angle.
     *  \details The tag of an angle should not be considered stable. Use with
     *  caution.
     *  \param tag the tag to set. */
    void SetTag(int32_t tag);
    
    /*! \brief Get the atoms of the angle.
     *  \return triple of the atoms of the angle. */
    const AngleAtoms& GetAtoms() const;
    
    /*! \brief Number of atoms this angle contains.
     *  \return 3. */
    constexpr int64_t NumAtoms() const { return  3; }
    
    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the angle in the container of angles
     *  of the molecule it is a part of. If the angle is not a real angle, eg
     *  empty data variable, the returned value is -1.
     *  \return the index of the angle. */
    int64_t GetIndex() const;
    
    /*! \brief Get the assigned type.
     *  \return the assigned ff type for the angle. */
    const FFAngle& GetType() const;
    
    /*! \brief Assign an FFAngle type.
     *  \param type the FFAngle type to assign. */
    void SetType(const FFAngle& type);
    
    bool HasType() const;
    
  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
    
    //! The molecule this angle is assigned to.
    wMolecule _mol;
    //! Tag (unstable)
    uint32_t _tag;
    //! MM angle type
    FFAngle _type;
    //! \brief Atoms which make up the angle.
    AngleAtoms _atms;
  };
  
  std::ostream& operator<<(std::ostream& os, const Angle& ang);
}

#endif  /* INDIGOX_CLASSES_ANGLE_HPP */
