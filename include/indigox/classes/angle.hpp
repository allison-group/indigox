/*! \file angle.hpp */
#include <array>
#include <memory>
#include <string>

#include "../utils/common.hpp"
#include "../utils/fwd_declares.hpp"
#include "../utils/counter.hpp"
#include "../utils/triple.hpp"

#ifndef INDIGOX_CLASSES_ANGLE_HPP
#define INDIGOX_CLASSES_ANGLE_HPP

namespace indigox {
  class Angle
  : public utils::IXCountableObject<Angle>,
  public std::enable_shared_from_this<Angle> {
    //! \brief Friendship allows Molecule to create new angles.
    friend class indigox::Molecule;
    //! \brief Friendship allows Angle to be tested.
    friend struct indigox::test::TestAngle;
    //! \brief Friendship allows serialisation.
    friend class cereal::access;
    
  private:
    //! \brief Container for storing Atom reference assigned to an Angle.
    using AngleAtoms = std::array<wAtom, 3>;
    
  public:
    // iterator aliases
    //! \brief Iterator over IXAtom references stored on an IXAngle.
    using AngleAtomIter = AngleAtoms::const_iterator;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<Angle>& construct,
                                   const uint32_t version);
    
    void Clear();
    
  public:
    Angle() = delete;  // no default constructor
    
    /*! \brief Normal constructor.
     *  \details Creates an angle between three atoms, linking it to the given
     *  Molecule. Atom \p b is always the central atom of the angle.
     *  \param a,b,c the three atoms to construct an angle between.
     *  \param m the molecule to assign the angle to. */
    Angle(Atom& a, Atom& b, Atom& c, Molecule& m);
    
    //! \brief Destructor
    ~Angle() { }
    
    /*! \brief Tag of the angle.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the angle, use the IXAngle::GetUniqueID().
     *  \return the tag assigned to the angle. */
    inline uint32_t GetTag() const { return _tag; }
    
    /*! \brief Molecule this angle is associated with.
     *  \details The returned shared_ptr is empty if the angle is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this angle. */
    Molecule& GetMolecule() const { return *_mol.lock(); }
    
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
    inline void SetTag(uint32_t tag) { _tag = tag; }
    
    /*! \brief Switch the begin and end atoms of the angle.
     *  \details The beginning atom will become the end atom and vice versa. The
     *  central atom will always remain the central atom. */
    inline void SwapOrder() { std::swap(_atms[0], _atms[2]); }
    
    /*! \brief Get the atoms of the angle.
     *  \return triple of the atoms of the angle. */
    inline stdx::triple<Atom&, Atom&, Atom&> GetAtoms() const {
      return {*_atms[0].lock(), *_atms[1].lock(), *_atms[2].lock()};
    }
    
    /*! \brief Number of atoms this angle contains.
     *  \return 3. */
    size_t NumAtoms() const { return  _atms.size(); }
    
    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the angle in the container of angles
     *  of the molecule it is a part of. If the molecule is dead, the index
     *  returned is the tag of the angle.
     *  \return the index of the angle. */
    size_t GetIndex() const;
    
    /*! \brief Get the assigned type.
     *  \return the assigned ff type for the angle. */
    FFAngle& GetType() const { return *_type.lock(); }
    
    /*! \brief Assign an FFAngle type.
     *  \param type the FFAngle type to assign. */
    void SetType(FFAngle& type);
    
    bool HasType() const { return !_type.expired(); }
    
  private:
    //! The molecule this angle is assigned to.
    wMolecule _mol;
    //! Tag (unstable)
    uint32_t _tag;
    //! MM angle type
    wFFAngle _type;
    //! \brief Atoms which make up the angle.
    AngleAtoms _atms;
  };
  
  std::ostream& operator<<(std::ostream& os, const Angle& ang);
}

#endif  /* INDIGOX_CLASSES_ANGLE_HPP */
