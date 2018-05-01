/*! \file atom.hpp */

#ifndef INDIGOX_CLASSES_ATOM_HPP
#define INDIGOX_CLASSES_ATOM_HPP

#include <memory>
#include <string>
#include <vector>

#include "../utils/counter.hpp"
#include "../utils/numerics.hpp"

namespace indigox {
  class IXAtom;
  class IXBond;
  class IXAngle;
  class IXDihedral;
  class IXMolecule;
  
  class IXElement;
  
  //! \brief shared_ptr for normal use of the IXAtom class.
  typedef std::shared_ptr<IXAtom> Atom;
  typedef std::shared_ptr<IXBond> Bond;
  typedef std::shared_ptr<IXAngle> Angle;
  typedef std::shared_ptr<IXDihedral> Dihedral;
  typedef std::shared_ptr<IXMolecule> Molecule;
  typedef std::shared_ptr<IXElement> Element;
  
  /*! \brief weak_ptr for non-ownership reference to the IXAtom class.
   *  \details Intended for internal use only. */
  typedef std::weak_ptr<IXAtom> _Atom;
  typedef std::weak_ptr<IXBond> _Bond;
  typedef std::weak_ptr<IXAngle> _Angle;
  typedef std::weak_ptr<IXDihedral> _Dihedral;
  typedef std::weak_ptr<IXMolecule> _Molecule;
  typedef std::weak_ptr<IXElement> _Element;
  
  //! \cond
  // Temporary defintion of Vec3 struct. Will make proper math stuff sometime.
  struct Vec3 {
    float_ x = 0.0, y = 0.0, z = 0.0;
  };
  //! \endcond
  
  class IXAtom
  : public utils::CountableObject<IXAtom>,
  public std::enable_shared_from_this<IXAtom> {
  
  private:
    // Typedefs
    //! \brief Container for storing IXBond references.
    typedef std::vector<_Bond> AtomBonds;
    //! \brief Container for storing IXAngle references.
    typedef std::vector<_Angle> AtomAngles;
    //! \brief Container for storing IXDihedral references.
    typedef std::vector<_Dihedral> AtomDihedrals;
  public:  // Make the iterator typedefs public for easier external usage
    //! \brief Iterator over IXBond references stored on an IXAtom
    typedef AtomBonds::iterator AtomBondIter;
    //! \brief Iterator over IXAngle references stored on an IXAtom
    typedef AtomAngles::iterator AtomAngleIter;
    //! \brief Iterator over IXDihedral references stored on an IXAtom
    typedef AtomDihedrals::iterator AtomDihedralIter;
    
  public:
    //! \brief Enum for the different types of atom stereochemistry
    enum class Stereo {
      UNDEFINED, //!< No defined stereochemistry.
      ACHIRAL,   //!< Defined as no stereochemistry.
      R,         //!< Has R stereochemistry.
      S,         //!< Has S stereochemistry.
    };

    /*! \brief Default constructor.
     *  \details Though allowed, it is not recommended to construct IXAtom
     *  instances directly. Rather, do so through the IXMolecule::NewAtom
     *  methods. */
    IXAtom();
    
    /*! \brief Normal constructor.
     *  \details Links the constructed atom to the given Molecule, though no
     *  bookkeeping is performed so the molecule does not know about it. Though
     *  allowed, it is not recommended to construct IXAtom instances directly.
     *  Rather, do so through the IXMolecule::NewAtom methods.
     *  \param m the molecule to assign this atom to. */
    IXAtom(Molecule m);
    
    //! \cond
    ~IXAtom() = default;
    //! \endcond
    
    /*! \brief Element of the atom.
     *  \return the element of this atom. */
    Element GetElement() const;
    
    /*! \brief Formal charge on the atom.
     *  \return the formal charge on the atom. */
    int_ GetFormalCharge() const { return _fc; }
    
    /*! \brief Partial atomic charge on the atom.
     *  \return the partial atomic charge. */
    float_ GetPartialCharge() const { return _partial; }
    
    /*! \brief Tag of the atom.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the atom, use IXAtom::GetUniqueID.
     *  \return the tag assigned to the atom. */
    uint_ GetTag() const { return _tag; };
    
    /*! \brief Get number of implicit hydrogens.
     *  \return the number of implicit hydrogens in the atom. */
    uint_ GetImplicitCount() const { return _implicitH; }
    
    /*! \brief Molecule this atom is associated with.
     *  \return the molecule associated with this atom.
     *  \throw std::logic_error Error if the atom was never assigned to a
     *  molecule, or if the assigned molecule has been deleted. */
    Molecule GetMolecule() const;
    
    /*! \brief Atom name.
     *  \return name of the atom. */
    string_ GetName() const { return _name; }
    
    /*! \brief Atom x position.
     *  \return the x coordinate of this atom. */
    float_ GetX() const { return _pos.x; }
    
    /*! \brief Atom y position.
     *  \return the y coordinate of this atom. */
    float_ GetY() const { return _pos.y; }
    
    /*! \brief Atom z position.
     *  \return the z coordinate of this atom. */
    float_ GetZ() const { return _pos.z; }
    
    /*! \brief Vector of the atom's position.
     *  \return the atoms position. */
    Vec3 GetVector() const { return _pos; }
    
    /*! \brief String representation of the atom.
     *  \details The returned string is of the form: Atom(NAME, SYMBOL).
     *  \return a string representation of the atom. */
    string_ ToString();
    
    /*! \brief Set the element of this atom.
     *  \param e the element to set to. */
    void SetElement(Element e) { _elem = e; }
    
    /*! \brief Set the element of this atom.
     *  \param e the name or atomic symbol of the element to set. */
    void SetElement(string_ e);
    
    /*! \brief Set the element of this atom.
     *  \param e the atomic number of the element to set. */
    void SetElement(uint_ e);
    
    /*! \brief Set the formal charge of this atom.
     *  \param q the formal charge value to set. */
    void SetFormalCharge(int_ q) { _fc = q; }
    
    /*! \brief Set the partial charge of this atom.
     *  \param q the partial charge value to set. */
    void SetPartialCharge(float_ q) { _partial = q; }
    
    /*! \brief Set the number of implicit hydrogens.
     *  \param h the number of implicit hydrogens to set. */
    void SetImplicitCount(uint_ h) { _implicitH = h; }
    
    /*! \brief Set the tag of this atom.
     *  \details The tag of an atom should not be considered stable. Use with
     *  caution.
     *  \param i the tag to set. */
    void SetTag(uint_ i) { _tag = i; }
    
    /*! \brief Set the molecule this atom is part of.
     *  \details No bookkeeping is performed, meaning the molecule is not
     *  informed that it now contains another atom. As such, this method is
     *  only intended for internal use.
     *  \param m the molecule to set. */
    void SetMolecule(Molecule m) { _mol = m; }
    
    /*! \brief Set the atom name.
     *  \param n name to set. */
    void SetName(string_ n) { _name = n; }
    
    /*! \brief Set the x position.
     *  \param x position to set. */
    void SetX(float_ x) { _pos.x = x; }
    
    /*! \brief Set the y position.
     *  \param y position to set. */
    void SetY(float_ y) { _pos.y = y; }
    
    /*! \brief Set the z position.
     *  \param z position to set. */
    void SetZ(float_ z) { _pos.z = z; }
    
    /*! \brief Set the x, y and z positions.
     *  \param x,y,z position to set. */
    void SetPosition(float_ x, float_ y, float_ z) {
      _pos.x = x; _pos.y = y; _pos.z = z;
    }
    
    /*! \brief Set the stereochemistry of an atomic center.
     *  \param s the stereochemistry to set. */
    void SetStereochemistry(Stereo s) { _stereo = s; }
    
    /*! \brief Set the aromaticity of an atom.
     *  \param a if the atom is aromatic or not. */
    void SetAromaticity(bool a) { _aromatic = a; }
    
    /*! \brief Get the stereochemistry of the atom.
     *  \return the stereochemistry of the atom. */
    Stereo GetStereochemistry() { return _stereo; }
    
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
    
    /*! \brief Number of bonds this atom is part of.
     *  \returns the number of assigned bonds. */
    size_ NumBonds() const { return _bonds.size(); }
    
    /*! \brief Number of angles this atom is a part of.
     *  \returns the number of assigned angles. */
    size_ NumAngles() const { return _angles.size(); }
    
    /*! \brief Number of dihedrals this atom is a part of.
     *  \returns the number of assigned dihedrals. */
    size_ NumDihedrals() const { return _dihedrals.size(); }
    
  private:
    //! The molecule this atom is assigned to.
    _Molecule _mol;
    //! The atoms element.
    _Element _elem;
    //! Formal charge.
    int_ _fc;
    //! Tag (unstable).
    uint_ _tag;
    //! Number of implicit hydrogens.
    uint_ _implicitH;
    //! Atoms name.
    string_ _name;
    //! Position vector.
    Vec3 _pos;
    //! Partial atomic charge.
    float_ _partial;
    //! Stereochemistry
    Stereo _stereo;
    //! Aromaticity
    bool _aromatic;
    
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
