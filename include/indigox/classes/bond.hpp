/*! \file bond.hpp */
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "../utils/counter.hpp"
#include "../utils/numerics.hpp"

#ifndef INDIGOX_CLASSES_BOND_HPP
#define INDIGOX_CLASSES_BOND_HPP

namespace indigox {
  class IXAtom;
  class IXBond;
  class IXAngle;
  class IXDihedral;
  class IXMolecule;
  
  typedef std::shared_ptr<IXAtom> Atom;
  //! \brief shared_ptr for normal use of the IXBond class.
  typedef std::shared_ptr<IXBond> Bond;
  typedef std::shared_ptr<IXAngle> Angle;
  typedef std::shared_ptr<IXDihedral> Dihedral;
  typedef std::shared_ptr<IXMolecule> Molecule;
  
  typedef std::weak_ptr<IXAtom> _Atom;
  /*! \brief weak_ptr for non-ownership reference to the IXBond class.
   *  \details Intended for internal use only. */
  typedef std::weak_ptr<IXBond> _Bond;
  typedef std::weak_ptr<IXAngle> _Angle;
  typedef std::weak_ptr<IXDihedral> _Dihedral;
  typedef std::weak_ptr<IXMolecule> _Molecule;
  
  class IXBond
  : public utils::CountableObject<IXBond>,
  public std::enable_shared_from_this<IXBond> {
  private:
    //! \brief Container for storing IXAtom references assigned to an IXBond.
    typedef std::array<_Atom, 2> BondAtoms;
    //! \brief Container for storing IXAngle references on an IXBond.
    typedef std::vector<_Angle> BondAngles;
    //! \brief Container for storing IXDihedral references on an IXBond.
    typedef std::vector<_Dihedral> BondDihedrals;
    
  public: // public iterator typedefs
    //! \brief Iterator over IXAtom references stored on an IXBond.
    typedef BondAtoms::iterator BondAtomIter;
    //! \brief Iterator over IXAngle references stored on an IXBond.
    typedef BondAngles::iterator BondAngleIter;
    //! \brief Iterator over IXDihedral references stored on an IXBond.
    typedef BondDihedrals::iterator BondDihedralIter;
    
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
    
  public:
    /*! \brief Default constructor.
     *  \details Though allowed, it is not recommended to construct IXBond
     *  instances directly. Rather, do so through the IXMolecule::NewBond
     *  methods. */
    IXBond();
    
    /*! \brief Normal constructor.
     *  \details Creates a bond between the two atoms, though no bookkeeping is
     *  performed, nor any sanity checks on the provided atoms. Tough allowed,
     *  it is not recommended to construct IXBond instances directly. Rather, do
     *  so through the IXMolecule::NewBond methods.
     *  \param a, b the atoms to construct a bonds between. */
    IXBond(Atom a, Atom b);
    
    //! \cond
    ~IXBond() = default;
    //! \endcond
    
    /*! \brief Tag of the bond.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the bond, use IXBond::GetUniqeID().
     *  \return the tag assigned to the bond. */
    uid_ GetTag() const { return _tag; }
    
    /*! \brief Molecule this bond is associated with.
     *  \details The returned shared_ptr is empty if the bond is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this bond. */
    Molecule GetMolecule() const { return _mol.lock(); }
    
    /*! \brief Get the bond order of the bond.
     *  \return the bond order. */
    Order GetOrder() const { return _order; }
    
    /*! \brief Get the source atom of the bond.
     *  \details The returned shared_ptr is empty if the source atom was never
     *  set or is no longer valid.
     *  \return the source atom of the bond. */
    Atom GetSourceAtom() const { return _atoms[0].lock(); }
    
    /*! \brief Get the target atom of the bond.
     *  \details The returned shared_ptr is empty if the target atom was never
     *  set or is no longer valid.
     *  \return the target atom of the bond. */
    Atom GetTargetAtom() const { return _atoms[1].lock(); }
    
    /*! \brief Get the aromaticity of a bond.
     *  \return if the bond is aromatic or not. */
    bool GetAromaticity() const { return _aromatic; }
    
    /*! \brief Get the stereochemistry of the bond.
     *  \return the stereochemistry of the bond. */
    Stereo GetStereochemistry() const { return _stereo; }
    
    /*! \brief String representation of the bond.
     *  \details If either of the atoms of the bond are not provided, the
     *  returned string is Bond(MALFORMED), otherwise it is of the form:
     *  Bond(NAME_SOURCE, NAME_TARGET).
     *  \return a string representation of the bond. */
    string_ ToString() const;
    
    /*! \brief Set the tag of this bond.
     *  \details The tag of a bond should not be considered stable. Use with
     *  caution.
     *  \param tag the tag to set. */
    void SetTag(uid_ tag) { _tag = tag; }
    
    /*! \brief Set the molecule this bond is part of.
     *  \details No bookkeeping is performed, meaning the molecule is not
     *  informed that it now contains another bond. As such, this method is only
     *  intended for internal use.
     *  \param mol the molecule to set. */
    void SetMolecule(Molecule mol) { _mol = mol; }
    
    /*! \brief Set the bond order.
     *  \param order the bond order to set. */
    void SetOrder(Order order) { _order = order; }
    
    /*! \brief Set the source atom of the bond.
     *  \details No bookkeeping is performed, so the atom does not know that it
     *  has been added to a new bond. As such, this method is intended for
     *  internal use only.
     *  \param atom the source atom to set. */
    void SetSourceAtom(Atom atom) { _atoms[0] = atom; }
    
    /*! \brief Set the target atom of the bond.
     *  \details No bookkeeping is performed, so the atom does not know that it
     *  has been added to a new bond. As such, this method is intended for
     *  internal use only.
     *  \param atom the target atom to set. */
    void SetTargetAtom(Atom atom) { _atoms[1] = atom; }
    
    /*! \brief Set the aromaticity.
     *  \param aromatic of the atom is aromatic or not. */
    void SetAromaticity(bool aromatic) { _aromatic = aromatic; }
    
    /*! \brief Set the stereochemistry of the bond.
     *  \param stereo the stereochemistry to set. */
    void SetStereochemistry(Stereo stereo) { _stereo = stereo; }
    
    /*! \brief Clear all information.
     *  \details Erases all information stored on the bond, and resets
     *  everything back to a just created state. No bookkeeping is performed. */
    void Clear();
    
    /*! \brief Get iterator access to the atoms of the bond.
     *  \details Intended for internal use only as the bond does not own any
     *  of the atoms being iterated over.
     *  \return a pair of iterators for the beginning and end of the atoms. */
    std::pair<BondAtomIter, BondAtomIter> GetAtomIters() {
      return std::make_pair(_atoms.begin(), _atoms.end());
    }
    
    /*! \brief Get iterator access to the angles the bond is involved in.
     *  \details Intended for internal use only as the bond does not own any
     *  of the angles being iterated over.
     *  \return a pair of iterators for the beginning and end of the angles. */
    std::pair<BondAngleIter, BondAngleIter> GetAngleIters() {
      return std::make_pair(_angles.begin(), _angles.end());
    }
    
    /*! \brief Get iterator access to the dihedrals the bond is involved in.
     *  \details Intended for internal use only as the bond does not own any
     *  of the dihedrals being iterated over.
     *  \return a pair of iterators for the beginning and end of the dihedrals. */
    std::pair<BondDihedralIter, BondDihedralIter> GetDihedralIters() {
      return std::make_pair(_dihedrals.begin(), _dihedrals.end());
    }
    
    /*! \brief Add an angle to this bond.
     *  \details No bookkeeping is performed, meaning that the angle is not
     *  checked to ensure it actually contains the current bond. As such, this
     *  method is intended for internal use only.
     *  \param ang the angle to add. */
    void AddAngle(Angle ang) { _angles.emplace_back(ang); }
    
    /*! \brief Add a dihedral to this bond.
     *  \details No bookkeeping is performed, meaning that the dihedral is not
     *  checked to ensure it actually contains the current bond. As such, this
     *  method is intended for internal use only.
     *  \param dihed the dihedral to add. */
    void AddDihedral(Dihedral dihed) { _dihedrals.emplace_back(dihed); }
    
    /*! \brief Number of atoms this bond is made of.
     *  \returns 2. */
    size_ NumAtoms() const { return _atoms.size(); }
    
    /*! \brief Number of angles this bond is a part of.
     *  \returns the number of assigned angles. */
    size_ NumAngles() const { return _angles.size(); }
    
    /*! \brief Number of dihedrals this bond is a part of.
     *  \returns the number of assigned dihedrals. */
    size_ NumDihedrals() const { return _dihedrals.size(); }
    
    /*! \brief Remove an angle from this bond.
     *  \details No bookkeeping is performed, meaning that the angle is not
     *  checked to ensure it actually contains the current bond. As such this
     *  method is intened for internal use only. If the angle is not part of
     *  this bond, nothing happens.
     *  \param ang the angle to remove. */
    void RemoveAngle(Angle ang);
    
    /*! \brief Remove a dihedral from this bond.
     *  \details No bookkeeping is performed, meaning that the dihedral is not
     *  checked to ensure it actually contains the current bond. As such this
     *  method is intened for internal use only. If the dihedral is not part of
     *  this bond, nothing happens.
     *  \param dihed the dihedral to remove. */
    void RemoveDihedral(Dihedral dihed);
    
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
    BondAtoms _atoms;
    //! \brief Angles which the bond is part of.
    BondAngles _angles;
    //! \brief Dihedrals which the bond is part of.
    BondDihedrals _dihedrals;
  };
  
}

#endif /* INDIGOX_CLASSES_BOND_HPP */
