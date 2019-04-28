/*! \file atom.hpp */

#ifndef INDIGOX_CLASSES_ATOM_HPP
#define INDIGOX_CLASSES_ATOM_HPP

#include "../utils/atomic_coordinates.hpp"
#include "../utils/fwd_declares.hpp"

#include <memory>
#include <vector>

namespace indigox {
  class Atom {
    //! \brief Friendship allows IXMolecule to create new atoms.
    friend class indigox::Molecule;
    //! \brief Friendship allows IXAtom to be serialised.
    friend class cereal::access;

  public:
    //! \brief Container for storing Bond references.
    using AtomBonds = std::vector<Bond>;
    //! \brief Container for storing IXAngle references.
    using AtomAngles = std::vector<Angle>;
    //! \brief Container for storing IXDihedral references.
    using AtomDihedrals = std::vector<Dihedral>;

  public:
    //! \brief Enum for the different types of atom stereochemistry
    enum class Stereo {
      UNDEFINED, //!< No defined stereochemistry.
      ACHIRAL,   //!< Defined as no stereochemistry.
      R,         //!< Has R stereochemistry.
      S,         //!< Has S stereochemistry.
    };

  private:
    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Atom);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Atom, atm);

  private:
    /*! \brief Normal constructor.
     *  \details Links the constructed atom to the given Molecule.
     *  \param m the molecule to assign this atom to. */
    Atom(const Molecule &molecule, const Element &element, std::string name);

  public:
    /*! \brief Number of valid bonds this atom is part of.
     *  \returns the number of valid assigned bonds. */
    int64_t NumBonds() const;

    /*! \brief Number of valid angles this atom is a part of.
     *  \returns the number of valid assigned angles. */
    int64_t NumAngles() const;

    /*! \brief Number of valid dihedrals this atom is a part of.
     *  \returns the number of valid assigned dihedrals. */
    int64_t NumDihedrals() const;

    bool HasType() const;

    /*! \brief Element of the atom.
     *  \return the element of this atom. */
    const Element &GetElement() const;

    /*! \brief Formal charge on the atom.
     *  \return the formal charge on the atom. */
    int32_t GetFormalCharge() const;

    /*! \brief Partial atomic charge on the atom.
     *  \return the partial atomic charge. */
    double GetPartialCharge() const;

    /*! \brief Tag of the atom.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the atom, use IXAtom::GetUniqueID.
     *  \return the tag assigned to the atom. */
    int64_t GetTag() const;

    int64_t GetID() const;

    int32_t GetResidueID();

    std::string GetResidueName();

    int32_t GetChargeGroupID() const;

    /*! \brief Get number of implicit hydrogens.
     *  \return the number of implicit hydrogens in the atom. */
    int32_t GetImplicitCount() const;

    int32_t NumHydrogenBonds() const;
    int32_t NumHeavyAtomBonds() const;

    /*! \brief Molecule this atom is associated with.
     *  \details The returned shared_ptr is empty of the atom is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this atom. */
    const Molecule &GetMolecule() const;

    /*! \brief Atom name.
     *  \return name of the atom. */
    const std::string &GetName() const;

    /*! \brief Atom x position.
     *  \return the x coordinate of this atom. */
    double GetX() const;

    /*! \brief Atom y position.
     *  \return the y coordinate of this atom. */
    double GetY() const;

    /*! \brief Atom z position.
     *  \return the z coordinate of this atom. */
    double GetZ() const;

    /*! \brief Get the stereochemistry of the atom.
     *  \return the stereochemistry of the atom. */
    Stereo GetStereochemistry() const;

    /*! \brief Vector of the atom's position.
     *  \return the atoms position. */
    Coordinates GetPosition() const;

    /*! \brief Add an implicit hydrogen.
     *  \return the new number of implicit hydrogens in the atom. */
    int32_t AddImplicitHydrogen();

    /*! \brief Remove an implicit hydrogen.
     *  \details If GetImplicitCount() == 0, no hydrogen is removed.
     *  \return the new number of implicit hydrogens in the atom. */
    int32_t RemoveImplicitHydrogen();

    /*! \brief Set the element of this atom.
     *  \param e the element to set to. */
    void SetElement(const Element &e);

    /*! \brief Set the element of this atom.
     *  \param e the name or atomic symbol of the element to set. */
    void SetElement(std::string e);

    /*! \brief Set the element of this atom.
     *  \param e the atomic number of the element to set. */
    void SetElement(int32_t e);

    /*! \brief Set the formal charge of this atom.
     *  \param q the formal charge value to set. */
    void SetFormalCharge(int32_t q);

    void SetChargeGroupID(int32_t id);

    /*! \brief Set the partial charge of this atom.
     *  \param q the partial charge value to set. */
    void SetPartialCharge(double q);

    /*! \brief Set the number of implicit hydrogens.
     *  \param h the number of implicit hydrogens to set. */
    void SetImplicitCount(int32_t h);

    /*! \brief Set the tag of this atom.
     *  \details The tag of an atom should not be considered stable. Use with
     *  caution.
     *  \param i the tag to set. */
    void SetTag(int32_t i);

    /*! \brief Set the atom name.
     *  \param n name to set. */
    void SetName(std::string n);

    /*! \brief Set the x, y and z positions.
     *  \param x,y,z position to set. */
    void SetPosition(double x, double y, double z);

    /*! \brief Set the stereochemistry of an atomic center.
     *  \param s the stereochemistry to set. */
    void SetStereochemistry(Stereo s);

  private:
    /*! \brief Add a bond to this atom.
     *  \details Assumes that the bond is not already added to the atom.
     *  \param b the bond to add.
     *  \return if the bond was added or not. */
    void AddBond(const Bond &b);

    /*! \brief Add an angle to this atom.
     *  \details Assumes that the angle is not already added to the atom.
     *  \param a the angle to add.
     *  \return if the angle was added or not. */
    void AddAngle(const Angle &a);

    /*! \brief Add a dihedral to this atom.
     *  \details Assumes that the dihedral is not already added to the atom.
     *  \param d the dihedral to add.
     *  \return if the dihedral was added or not. */
    void AddDihedral(const Dihedral &d);

    /*! \brief Remove a bond from this atom.
     *  \details Assumes that the bond is added to the atom.
     *  \param b the bond to remove.
     *  \return if the bond was removed or not. */
    void RemoveBond(const Bond &b);

    /*! \brief Remove an angle from this atom.
     *  \details Assumes that the angle is added to the atom.
     *  \param a the angle to remove.
     *  \return if the angle was removed or not. */
    void RemoveAngle(const Angle &a);

    /*! \brief Remove a dihedral from this atom.
     *  \details Assumes that the dihedral is added to the atom.
     *  \param d the dihedral to remove.
     *  \return if the dihedral was removed or not. */
    void RemoveDihedral(const Dihedral &d);

  public:
    const AtomBonds &GetBonds() const;
    const AtomAngles &GetAngles() const;
    const AtomDihedrals &GetDihedrals() const;

    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the atom in the container of atoms
     *  of the molecule it is a part of. If the molecule is dead, the index
     *  returned is the tag of the atom.
     *  \return the index of the atom. */
    int64_t GetIndex() const;

    /*! \brief Get the FF type of the atom.
     *  \return the force field type of the atom. */
    const FFAtom &GetType() const;

    /*! \brief Set the FF type of the atom.
     *  \param type the type to set. */
    void SetType(const FFAtom &type);

  private:
    void Reset();

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
  };

  //! \brief Type for the stereochemistry enum of an atom.
  using AtomStereo = indigox::Atom::Stereo;
} // namespace indigox

#endif /* INDIGOX_CLASSES_ATOM_HPP */
