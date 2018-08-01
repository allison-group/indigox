/*! \file bond.hpp */
#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../utils/common.hpp"
#include "../utils/counter.hpp"
#include "../utils/numerics.hpp"

#ifndef INDIGOX_CLASSES_BOND_HPP
#define INDIGOX_CLASSES_BOND_HPP

namespace indigox {
  class IXAtom;
  class IXBond;
  class IXMolecule;
  namespace test { struct TestBond; }
  
  using Atom = std::shared_ptr<IXAtom>;
  //! \brief shared_ptr for normal use of the IXBond class.
  using Bond = std::shared_ptr<IXBond>;
  using Molecule = std::shared_ptr<IXMolecule>;
  
  using _Atom = std::weak_ptr<IXAtom>;
  /*! \brief weak_ptr for non-ownership reference to the IXBond class.
   *  \details Intended for internal use only. */
  using _Bond = std::weak_ptr<IXBond>;
  using _Molecule = std::weak_ptr<IXMolecule>;
  
  class IXBond
  : public utils::IXCountableObject<IXBond>,
  public std::enable_shared_from_this<IXBond> {
    //! \brief Friendship allows IXMolecule to create new bonds.
    friend class indigox::IXMolecule;
    //! \brief Friendship allows IXBond to be tested.
    friend struct indigox::test::TestBond;
    //! \brief Friendship allows serialisation.
    friend class cereal::access;
    
  private:
    //! \brief Container for storing IXAtom references assigned to an IXBond.
    using BondAtoms = std::array<_Atom, 2>;
    
  public: // public iterator aliases
    //! \brief Iterator over IXAtom references stored on an IXBond.
    using BondAtomIter = BondAtoms::const_iterator;
    
  public:
    //! \brief Enum for the different possible bond stereochemistry states.
    enum class Stereo {
      UNDEFINED,  //!< No defined stereochemistry.
      NONE,       //!< Defined as no stereochemistry.
      E,          //!< E isomer.
      Z,          //!< Z isomer.
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
    /*! \brief Normal constructor.
     *  \details Creates a bond between the two atoms, linking it to the given
     *  Molecule.
     *  \param a, b the atoms to construct a bonds between.
     *  \param m the molecule to assign the bond to. */
    IXBond(Atom a, Atom b, Molecule m);
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXBond>& construct,
                                   const uint32_t version);
    
  public:
    IXBond() = delete;  // no default constructor
    
    //! \brief Destructor
    ~IXBond() { }
    
    /*! \brief Tag of the bond.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the bond, use IXBond::GetUniqeID().
     *  \return the tag assigned to the bond. */
    inline uid_ GetTag() const { return _tag; }
    
    /*! \brief Molecule this bond is associated with.
     *  \details The returned shared_ptr is empty if the bond is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this bond. */
    inline Molecule GetMolecule() const { return _mol.lock(); }
    
    /*! \brief Get the bond order of the bond.
     *  \return the bond order.
     *  \todo Obtain order based on molecule data. */
    inline Order GetOrder() const { return _order; }
    
    /*! \brief Get the source atom of the bond.
     *  \details The returned shared_ptr is empty if the source atom is no
     *  longer valid.
     *  \return the source atom of the bond. */
    inline Atom GetSourceAtom() const { return _atms[0].lock(); }
    
    /*! \brief Get the target atom of the bond.
     *  \details The returned shared_ptr is empty if the target atom is no
     *  longer valid.
     *  \return the target atom of the bond. */
    inline Atom GetTargetAtom() const { return _atms[1].lock(); }
    
    /*! \brief Get the aromaticity of a bond.
     *  \return if the bond is aromatic or not.
     *  \todo Obtain aromaticity based on molecule data. */
    inline bool GetAromaticity() const { return _aromatic; }
    
    /*! \brief Get the stereochemistry of the bond.
     *  \return the stereochemistry of the bond.
     *  \todo Obtain stereochemistry based on molecule data. */
    inline Stereo GetStereochemistry() const { return _stereo; }
    
    /*! \brief String representation of the bond.
     *  \details If either of the atoms of the bond are not valid, the
     *  returned string is Bond(MALFORMED), otherwise it is of the form:
     *  Bond(NAME_SOURCE, NAME_TARGET).
     *  \return a string representation of the bond. */
    string_ ToString() const;
    
    /*! \brief Set the tag of this bond.
     *  \details The tag of a bond should not be considered stable. Use with
     *  caution.
     *  \param tag the tag to set. */
    inline void SetTag(uid_ tag) { _tag = tag; }
    
    /*! \brief Set the bond order.
     *  \param order the bond order to set. */
    inline void SetOrder(Order order) { _order = order; }
    
    /*! \brief Switch the two atoms of the bond.
     *  \details the source atom will become the target atom and vice versa. */
    inline void SwapSourceTarget() { std::swap(_atms[0], _atms[1]); }
    
    /*! \brief Set the aromaticity.
     *  \param aromatic of the atom is aromatic or not. */
    inline void SetAromaticity(bool aromatic) { _aromatic = aromatic; }
    
    /*! \brief Set the stereochemistry of the bond.
     *  \param stereo the stereochemistry to set. */
    inline void SetStereochemistry(Stereo stereo) { _stereo = stereo; }
    
    /*! \brief Get the atoms of the bond.
     *  \return pair of the atoms of the bond. */
    inline std::pair<Atom, Atom> GetAtoms() const {
      return std::make_pair(_atms[0].lock(), _atms[1].lock());
    }
    
  private:
    /*! \brief Clear all information.
     *  \details Erases all information stored on the bond, and resets
     *  everything back to a just created state. */
    void Clear();
    
  public:
    /*! \brief Number of atoms this bond is between.
     *  \returns 2. */
    size_ NumAtoms() const { return _atms.size(); }
    
    size_ GetIndex() const;
    
  private:
    //! The molecule this bond is assigned to.
    _Molecule _mol;
    //! Tag (unstable).
    uid_ _tag;
    //! Bond order.
    Order _order;
    //! Aromaticity.
    bool _aromatic;
    //! Stereochemistry.
    Stereo _stereo;
    
    //! \brief Atoms which make up the bond.
    BondAtoms _atms;
  };
  
  /*! \brief Print a Bond to an output stream.
   *  \details The printed string is of the form: Bond(SOURCE, TARGET).
   *  \param os the output stream to print to.
   *  \param bnd the Bond to print.
   *  \return the output stream after printing. */
  std::ostream& operator<<(std::ostream& os, const IXBond& bnd);
  inline std::ostream& operator<<(std::ostream& os, const Bond& bnd) {
    return bnd ? os << *bnd : os;
  }
  
  //! Type for the stereochemistry enum of a bond.
  using BondStereo = indigox::IXBond::Stereo;
  //! Type for the order of a bond.
  using BondOrder = indigox::IXBond::Order;
  
}

#endif /* INDIGOX_CLASSES_BOND_HPP */
