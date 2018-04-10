/** @file atom.hpp
 *  @brief Atom declaration
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#ifndef INDIGOX_CLASSES_ATOM_HPP
#define INDIGOX_CLASSES_ATOM_HPP

#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "periodictable.hpp"
#include "../utils/counter.hpp"

namespace indigox {

  // Related typedefs
  //! \cond
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
  
  typedef std::vector<_Bond> AtomBonds;
  typedef std::vector<_Angle> AtomAngles;
  typedef std::vector<_Dihedral> AtomDihedrals;
  typedef AtomBonds::iterator AtomBondIter;
  typedef AtomAngles::iterator AtomAngleIter;
  typedef AtomDihedrals::iterator AtomDihedralIter;
  
  enum ATOMSTEREO {
    ACHIRAL,
    CHIRAL_R,
    CHIRAL_S
  };
  
  // Temporary defintion of Vec3 struct. Will make proper math stuff sometime.
  struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
  };
  //! \endcond
  
  class IXAtom
  : public utils::CountableObject<IXAtom>,
  public std::enable_shared_from_this<IXAtom> {
  
  public:

    /*! \brief Default constructor.
     *  \details Though allowed, it is not recommended to construct IXAtom
     *  instances directly. Rather, do so through the IXMolecule::NewAtom
     *  methods. */
    IXAtom();
    
    /*! \brief Normal constructor, links IXAtom to a IXMolecule.
     *  \details Though allowed, it is not recommended to construct IXAtom
     *  instances directly. Rather, do so through the IXMolecule::NewAtom
     *  methods.
     *  \param m the molecule to assign this atom to. */
    IXAtom(Molecule m);
    
    //! \brief Destructor
    ~IXAtom() = default;
    
    /*! \brief Element of the atom.
     *  \return the element of this atom. */
    Element GetElement() const {
      if(_elem.expired())
        return IXPeriodicTable::GetInstance()->GetUndefinedElement();
      return _elem.lock();
    }
    
    /*! \brief Formal charge on the atom.
     *  \return the formal charge on the atom. */
    int GetFormalCharge() const { return _fc; }
    
    /*! \brief Partial atomic charge on the atom.
     *  \return the partial atomic charge. */
    double GetPartialCharge() const { return _partial; }
    
    /*! \brief Index of the atom.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the atom, use IXAtom::GetUniqueID.
     *  \return the index assigned to the atom. */
    unsigned int GetIndex() const { return _idx; };
    
    /*! \brief Get number of implicit hydrogens.
     *  \return the number of implicit hydrogens in the atom. */
    unsigned int GetImplicitCount() const { return _implicitH; }
    
    /*! \brief Molecule this atom is associated with.
     *  \return the molecule associated with this atom.
     *  \throw std::logic_error Error if the atom was never assigned to a
     *  molecule, or if the assigned molecule has been deleted. */
    Molecule GetMolecule() const {
      if(_mol.expired())
        throw std::logic_error("Atom not assigned to a valid molecule");
      return _mol.lock();
    }
    
    /*! \brief Atom name.
     *  \return name of the atom. */
    std::string GetName() const { return _name; }
    
    /*! \brief Atom x position.
     *  \return the x coordinate of this atom. */
    double GetX() const { return _pos.x; }
    
    /*! \brief Atom y position.
     *  \return the y coordinate of this atom. */
    double GetY() const { return _pos.y; }
    
    /*! \brief Atom z position.
     *  \return the z coordinate of this atom. */
    double GetZ() const { return _pos.z; }
    
    /*! \brief Vector of the atom's position.
     *  \return the atoms position. */
    Vec3 GetVector() const { return _pos; }
    
    /*! \brief String representation of the atom.
     *  \details The representation returned depends on if indigox was compiled
     *  in debug mode or not. In release mode, the string is of the form:
     *  Atom(NAME, ELEMENT). In debug mode, the string is of the form:
     *  Atom(NAME-IDX, ELEMENT, X.X, Y.Y, Z.Z).
     *  \return a string representation of the atom. */
    std::string ToString();
    
    /*! \brief Set the element of this atom.
     *  \param e the element to set to. */
    void SetElement(Element e) { _elem = e; }
    
    /*! \brief Set the element of this atom.
     *  \param e the name or atomic symbol of the element to set. */
    void SetElement(std::string e) {
      _elem = IXPeriodicTable::GetInstance()->GetElement(e);
    }
    
    /*! \brief Set the element of this atom.
     *  \param e the atomic number of the element to set. */
    void SetElement(unsigned int e)  {
      _elem = IXPeriodicTable::GetInstance()->GetElement(e);
    }
    
    /*! \brief Set the formal charge of this atom.
     *  \param q the formal charge value to set. */
    void SetFormalCharge(int q) { _fc = q; }
    
    /*! \brief Set the partial charge of this atom.
     *  \param q the partial charge value to set. */
    void SetPartialCharge(double q) { _partial = q; }
    
    /*! \brief Set the number of implicit hydrogens.
     *  \param h the number of implicit hydrogens to set. */
    void SetImplicitCount(unsigned int h) { _implicitH = h; }
    
    /*! \brief Set the index of this atom.
     *  \details The index of an atom should not be considered stable. Use with
     *  caution.
     *  \param i the index to set. */
    void SetIndex(unsigned int i) { _idx = i; }
    
    /*! \brief Set the molecule this atom is part of.
     *  \details No bookkeeping is performed, meaning the molecule is not
     *  informed that it now contains another atom. As such, this method is
     *  only intended for internal use.
     *  \param m the molecule to set. */
    void SetMolecule(Molecule m) { _mol = m; }
    
    /*! \brief Set the atom name.
     *  \param n name to set. */
    void SetName(std::string n) { _name = n; }
    
    /*! \brief Set the x position.
     *  \param x position to set. */
    void SetX(double x) { _pos.x = x; }
    
    /*! \brief Set the y position.
     *  \param y position to set. */
    void SetY(double y) { _pos.y = y; }
    
    /*! \brief Set the z position.
     *  \param z position to set. */
    void SetZ(double z) { _pos.z = z; }
    
    /*! \brief Set the x, y and z positions.
     *  \param x,y,z position to set. */
    void SetPosition(double x, double y, double z) {
      _pos.x = x; _pos.y = y; _pos.z = z;
    }
    
    /*! \brief Set the stereochemistry of an atomic center.
     *  \param s the stereochemistry to set. */
    void SetStereochemistry(ATOMSTEREO s) { _stereo = s; }
    
    /*! \brief Set the aromaticity of an atom.
     *  \param a if the atom is aromatic or not. */
    void SetAromaticity(bool a) { _aromatic = a; }
    
    /*! \brief Get the stereochemistry of the atom.
     *  \return the stereochemistry of the atom. */
    ATOMSTEREO GetStereochemistry() { return _stereo; }
    
    /*! \brief Get the aromaticity of an atom.
     *  \return if the atom is aromatic or not. */
    bool GetAromaticity() { return _aromatic; }
    
    /*! \brief Add a bond to this atom.
     *  \details No bookkeeping is performed, meaning that the bond is not
     *  checked to ensure it actually contains the current atom. As such, this
     *  method is only intended for internal use.
     *  \param b the bond to add. */
    void AddBond(Bond b) { _bonds.emplace_back(b); }
    
    /*! \brief Add an angle to this atom.
     *  \details No bookkeeping is performed, meaning that the angle is not
     *  checked to ensure it actually contains the current atom. As such, this
     *  method is only intended for internal use.
     *  \param a the angle to add. */
    void AddAngle(Angle a) { _angles.emplace_back(a); }
    
    /*! \brief Add a dihedral to this atom.
     *  \details No bookkeeping is performed, meaning that the dihedral is not
     *  checked to ensure it actually contains the current atom. As such, this
     *  method is only intended for internal use.
     *  \param d the dihedral to add. */
    void AddDihedral(Dihedral d) { _dihedrals.emplace_back(d); }
    
    /*! \brief Remove a bond from this atom.
     *  \details No bookkeeping is performed, meaning that the bond is not
     *  checked to ensure it actually contains the current atom. As such, this
     *  method is only intended for internal use. If the bond is not part of
     *  this atom, nothing happens.
     *  \param b the bond to remove. */
    void RemoveBond(Bond b);
    
    /*! \brief Remove an angle from this atom.
     *  \details No bookkeeping is performed, meaning that the angle is not
     *  checked to ensure it actually contains the current atom. As such, this
     *  method is only intended for internal use. If the angle is not part of
     *  this atom, nothing happens.
     *  \param a the angle to remove. */
    void RemoveAngle(Angle a);
    
    /*! \brief Remove a dihedral from this atom.
     *  \details No bookkeeping is performed, meaning that the dihedral is not
     *  checked to ensure it actually contains the current atom. As such, this
     *  method is only intended for internal use. If the dihedral is not part of
     *  this atom, nothing happens.
     *  \param d the dihedral to remove. */
    void RemoveDihedral(Dihedral d);
    
    /*! \brief Clear all information.
     *  \details Erases all information stored on the atom. */
    void Clear();
    
    /*! \brief Get iterator access to the atom's bonds.
     *  \details Intended primarily for internal use as the iterators are to
     *  weak_ptrs.
     *  \returns a pair of iterators for the beginning and end of the bonds. */
    std::pair<AtomBondIter, AtomBondIter> GetBondIters() {
      return std::make_pair(_bonds.begin(), _bonds.end());
    }
    
    /*! \brief Get iterator access to the atom's angles.
     *  \details Intended primarily for internal use as the iterators are to
     *  weak_ptrs.
     *  \returns a pair of iterators for the beginning and end of the angles. */
    std::pair<AtomAngleIter, AtomAngleIter> GetAngleIters() {
      return std::make_pair(_angles.begin(), _angles.end());
    }
    
    /*! \brief Get iterator access to the atom's dihedrals.
     *  \details Intended primarily for internal use as the iterators are to
     *  weak_ptrs.
     *  \returns a pair of iterators for the beginning and end of the
     *  dihedrals. */
    std::pair<AtomDihedralIter, AtomDihedralIter> GetDihedralIters() {
      return std::make_pair(_dihedrals.begin(), _dihedrals.end());
    }
    
  private:
    //! The molecule this atom is assigned to.
    _Molecule _mol;
    //! The atoms element.
    _Element _elem;
    //! Formal charge.
    int _fc = 0;
    //! Index (unstable).
    unsigned int _idx = 0;
    //! Number of implicit hydrogens.
    unsigned int _implicitH = 0;
    //! Atoms name.
    std::string _name = "ATOM";
    //! Position vector.
    Vec3 _pos;
    //! Partial atomic charge.
    double _partial = 0.0;
    //! Stereochemistry
    ATOMSTEREO _stereo = ACHIRAL;
    //! Aromaticity
    bool _aromatic = false;
    
    //! MM type for atom
    // FFAtom _mmtype;
    
    //! Bonds the atom is part of
    AtomBonds _bonds;
    //! Angles the atom is part of
    AtomAngles _angles;
    //! Dihedrals the atom is part of
    AtomDihedrals _dihedrals;
    
  };
}

#endif /* INDIGOX_CLASSES_ATOM_HPP */
