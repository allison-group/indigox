/*! \file angle.hpp */

#include "../utils/fwd_declares.hpp"

#include <array>
#include <memory>

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
    void serialise(Archive &archive, const uint32_t version);

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Angle);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Angle, ang);

  private:
    /*! \brief Normal constructor.
     *  \details Creates an angle between three atoms, linking it to the given
     *  Molecule. Atom \p b is always the central atom of the angle.
     *  \param a,b,c the three atoms to construct an angle between.
     *  \param m the molecule to assign the angle to. */
    Angle(const Atom &a, const Atom &b, const Atom &c, const Molecule &m);

  public:
    /*! \brief Check if the angle has associated parameters.
     *  \return if there is a valid parameter. */
    bool HasType() const;

    /*! \brief Number of atoms this angle contains.
     *  \return 3. */
    constexpr int64_t NumAtoms() const { return 3; }

    /*! \brief Tag of the angle.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the angle, use the Angle::GetUniqueID()
     *  method.
     *  \return the tag assigned to the angle. */
    int64_t GetTag() const;

    int64_t GetID() const;

    /*! \brief Molecule this angle is associated with.
     *  \details The returned shared_ptr is empty if the angle is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this angle. */
    const Molecule &GetMolecule() const;

    /*! \brief Get the atoms of the angle.
     *  \return the atoms of the angle. */
    const AngleAtoms &GetAtoms() const;

    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the angle in the container of angles
     *  of the molecule it is a part of. If the angle is not a real angle, eg
     *  empty data variable, the returned value is -1.
     *  \return the index of the angle. */
    int64_t GetIndex() const;

    /*! \brief Get the assigned type.
     *  \return the assigned ff type for the angle. */
    const FFAngle &GetType() const;

    /*! \brief Assign an FFAngle type.
     *  \param type the FFAngle type to assign. */
    void SetType(const FFAngle &type);

  private:
    void Reset();

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
  };
} // namespace indigox

#endif /* INDIGOX_CLASSES_ANGLE_HPP */
