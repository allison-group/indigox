/** @file bond.hpp
 *  @brief Bond declaration
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */
#include <array>
#include <memory>
#include <string>

#include "../utils/counter.hpp"

#ifndef INDIGOX_CLASSES_BOND_HPP
#define INDIGOX_CLASSES_BOND_HPP

namespace indigox {
  
  class IXAtom;
  class IXBond;
  class IXAngle;
  class IXDihedral;
  class IXMolecule;
  typedef std::shared_ptr<IXAtom> Atom;
  typedef std::shared_ptr<IXBond> Bond;
  typedef std::shared_ptr<IXAngle> Angle;
  typedef std::shared_ptr<IXDihedral> Dihedral;
  typedef std::shared_ptr<IXMolecule> Molecule;
  typedef std::weak_ptr<IXAtom> _Atom;
  typedef std::weak_ptr<IXBond> _Bond;
  typedef std::weak_ptr<IXAngle> _Angle;
  typedef std::weak_ptr<IXDihedral> _Dihedral;
  typedef std::weak_ptr<IXMolecule> _Molecule;
  
  // Related typedefs
  typedef std::array<_Atom, 2> BondAtoms;
  typedef BondAtoms::iterator BondAtomIterator;
  typedef BondAtoms::const_iterator const_BondAtomIterator;
  
  enum BondOrder {
    SINGLE_BOND = 1,
    DOUBLE_BOND,
    TRIPLE_BOND,
    QUADRUPLE_BOND,
    AROMATIC_BOND,
    ONEANDAHALF_BOND,
    TWOANDAHALF_BOND,
    UNDEFINED_BOND
  };
  
  class IXBond
  : public utils::CountableObject<IXBond>,
  public std::enable_shared_from_this<IXBond> {
  
  public:
    
    /// @name Initialisation methods.
    
    /// @brief Default constructor.
    IXBond();
    
    /// @brief Normal use constructor.
    IXBond(Atom, Atom);
    
    /// @brief Destructor
    ~IXBond();
    
    /// @name Data retrival methods
    
    /// @returns the index of this bond.
    uid_t GetIndex() const;
    
    /// @returns the molecule this bond is part of.
    Molecule GetMolecule() const;
    
    /// @returns the bond order of this bond.
    BondOrder GetOrder() const;
    
    /// @returns the source atom of this bond.
    Atom GetSourceAtom() const;
    
    /// @returns the target atom of this bond.
    Atom GetTargetAtom() const;
    
    /// @returns a string representation of this bond.
    std::string ToString() const;
    
    /// @name Data modification methods
    
    /// @brief Sets the index of this bond.
    void SetIndex(uid_t);
    
    /// @brief Sets the molecule this bond is part of.
    void SetMolecule(Molecule);
    
    /// @brief Sets the bond order of this bond.
    void SetOrder(BondOrder);
    
    /// @brief Sets the source atom of this bond.
    void SetSourceAtom(Atom);
    
    /// @brief Sets the target atom of this bond.
    void SetTargetAtom(Atom);
    
    void Clear();
    
    /// @name Iterator methods
    
    /// @brief Sets the given iterator to the next position in atoms_.
    Atom Next(BondAtomIterator&);
    
    /// @brief Sets the given iterator to the start of atoms_.
    Atom Begin(BondAtomIterator&);
    
    /// @brief Sets the given iterator to the end of atoms_.
    BondAtomIterator BeginAtom();
    BondAtomIterator EndAtom();
    
  private:
    BondAtoms atoms_;
    std::weak_ptr<IXMolecule> mol_;
    BondOrder order_ = UNDEFINED_BOND;
    uid_t idx_;
  };
  
}

#endif /* INDIGOX_CLASSES_BOND_HPP */
