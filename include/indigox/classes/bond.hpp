/*! \file bond.hpp */
#include "../utils/fwd_declares.hpp"

#include <array>
#include <memory>

#ifndef INDIGOX_CLASSES_BOND_HPP
#define INDIGOX_CLASSES_BOND_HPP

namespace indigox {
  class Bond {
    //! \brief Friendship allows IXMolecule to create new bonds.
    friend class indigox::Molecule;
    //! \brief Friendship allows serialisation.
    friend class cereal::access;

  public:
    //! \brief Container for storing Atom references assigned to an Bond.
    using BondAtoms = std::array<Atom, 2>;

  public:
    //! \brief Enum for the different possible bond stereochemistry states.
    enum class Stereo {
      UNDEFINED, //!< No defined stereochemistry.
      NONE,      //!< Defined as no stereochemistry.
      E,         //!< E isomer.
      Z,         //!< Z isomer.
    };

    //! \brief Enum for the different possible bond orders
    enum class Order {
      UNDEFINED,   //!< No defined bond order.
      SINGLE,      //!< Bond order of 1.
      DOUBLE,      //!< Bond order of 2.
      TRIPLE,      //!< Bond order of 3.
      QUADRUPLE,   //!< Bond order of 4.
      AROMATIC,    //!< Aromatic bond.
      ONEANDAHALF, //!< Non-aromatic bond order of 1.5.
      TWOANDAHALF, //!< Bond order of 2.5
    };

  private:
    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Bond);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Bond, bnd);

  private:
    /*! \brief Normal constructor.
     *  \details Creates a bond between the two atoms, linking it to the given
     *  Molecule.
     *  \param a, b the atoms to construct a bonds between.
     *  \param m the molecule to assign the bond to. */
    Bond(const Atom &a, const Atom &b, const Molecule &m, Order o);

  public:
    /*! \brief Tag of the bond.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the bond, use IXBond::GetUniqeID().
     *  \return the tag assigned to the bond. */
    int64_t GetTag() const;
    int64_t GetID() const;

    /*! \brief Molecule this bond is associated with.
     *  \details The returned shared_ptr is empty if the bond is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this bond. */
    const Molecule &GetMolecule() const;

    /*! \brief Get the bond order of the bond.
     *  \return the bond order.
     *  \todo Obtain order based on molecule data. */
    Order GetOrder() const;

    /*! \brief Get the stereochemistry of the bond.
     *  \return the stereochemistry of the bond.
     *  \todo Obtain stereochemistry based on molecule data. */
    Stereo GetStereochemistry() const;

    /*! \brief Set the tag of this bond.
     *  \details The tag of a bond should not be considered stable. Use with
     *  caution.
     *  \param tag the tag to set. */
    void SetTag(int64_t tag);

    /*! \brief Set the bond order.
     *  \param order the bond order to set. */
    void SetOrder(Order order);

    /*! \brief Set the stereochemistry of the bond.
     *  \param stereo the stereochemistry to set. */
    void SetStereochemistry(Stereo stereo);

    /*! \brief Get the atoms of the bond.
     *  \return pair of the atoms of the bond. */
    const BondAtoms &GetAtoms() const;

  public:
    /*! \brief Number of atoms this bond is between.
     *  \returns 2. */
    constexpr int64_t NumAtoms() const { return 2; }

    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the bond in the container of bonds
     *  of the molecule it is a part of. If the molecule is dead, the index
     *  returned is the tag of the bond.
     *  \return the index of the bond. */
    int64_t GetIndex() const;

    /*! \brief Get the type of the bond.
     *  \return the type of the bond. */
    const FFBond &GetType() const;

    /*! \brief Set the type of the bond.
     *  \param type the type of bond to set. */
    void SetType(const FFBond &type);

    bool HasType() const;

    bool IsAmideBond() const;
    bool IsCarbonylBond() const;

  private:
    void Reset();

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
  };

  //! Type for the stereochemistry enum of a bond.
  using BondStereo = indigox::Bond::Stereo;
  //! Type for the order of a bond.
  using BondOrder = indigox::Bond::Order;

} // namespace indigox

#endif /* INDIGOX_CLASSES_BOND_HPP */
