/*! \file atom.hpp */

#ifndef INDIGOX_CLASSES_ATOM_HPP
#define INDIGOX_CLASSES_ATOM_HPP

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "../classes/periodictable.hpp"
#include "../utils/common.hpp"
#include "../utils/counter.hpp"
#include "../utils/numerics.hpp"

namespace indigox {
  class IXAtom;
  class IXBond;
  class IXAngle;
  class IXDihedral;
  class IXMolecule;
  class IXElement;
  class IXFFAtom;
  namespace test { struct TestAtom; }
  
  //! \brief shared_ptr for normal use of the IXAtom class.
  using Atom = std::shared_ptr<IXAtom>;
  using Bond = std::shared_ptr<IXBond>;
  using Angle = std::shared_ptr<IXAngle>;
  using Dihedral = std::shared_ptr<IXDihedral>;
  using Molecule = std::shared_ptr<IXMolecule>;
  using Element = std::shared_ptr<IXElement>;
  using FFAtom = std::shared_ptr<IXFFAtom>;
  
  /*! \brief weak_ptr for non-ownership reference to the IXAtom class.
   *  \details Intended for internal use only. */
  using _Atom = std::weak_ptr<IXAtom>;
  using _Bond = std::weak_ptr<IXBond>;
  using _Angle = std::weak_ptr<IXAngle>;
  using _Dihedral = std::weak_ptr<IXDihedral>;
  using _Molecule = std::weak_ptr<IXMolecule>;
  using _Element = std::weak_ptr<IXElement>;
  
  class IXAtom
  : public utils::IXCountableObject<IXAtom>,
  public std::enable_shared_from_this<IXAtom> {
    //! \brief Friendship allows IXMolecule to create new atoms.
    friend class indigox::IXMolecule;
    //! \brief Friendship allows IXAtom to be tested.
    friend struct indigox::test::TestAtom;
    //! \brief Friendship allows IXAtom to be serialised.
    friend class cereal::access;
    
  private:
    //! \brief Container for storing IXBond references.
    using AtomBonds = std::vector<_Bond>;
    //! \brief Container for storing IXAngle references.
    using AtomAngles = std::vector<_Angle>;
    //! \brief Container for storing IXDihedral references.
    using AtomDihedrals = std::vector<_Dihedral>;
  public:  // Make the iterator aliases public for easier external usage
    //! \brief Iterator over IXBond references stored on an IXAtom
    using AtomBondIter = AtomBonds::const_iterator;
    //! \brief Iterator over IXAngle references stored on an IXAtom
    using AtomAngleIter = AtomAngles::const_iterator;
    //! \brief Iterator over IXDihedral references stored on an IXAtom
    using AtomDihedralIter = AtomDihedrals::const_iterator;
    
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
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXAtom>& construct,
                                   const uint32_t version);
    
  public:
    IXAtom() = delete;  // default constructor for serialise access
    
    /*! \brief Normal constructor.
     *  \details Links the constructed atom to the given Molecule.
     *  \param m the molecule to assign this atom to. */
    IXAtom(Molecule m);
    
    //! \brief Destructor
    ~IXAtom() { };
    
    /*! \brief Element of the atom.
     *  \return the element of this atom. */
    inline Element GetElement() const {
      return _elem.expired() ? GetPeriodicTable()->GetUndefined() : _elem.lock();
    }
    
    /*! \brief Formal charge on the atom.
     *  \return the formal charge on the atom. */
    inline int_ GetFormalCharge() const { return _fc; }
    
    /*! \brief Partial atomic charge on the atom.
     *  \return the partial atomic charge. */
    inline float_ GetPartialCharge() const { return _partial; }
    
    /*! \brief Tag of the atom.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the atom, use IXAtom::GetUniqueID.
     *  \return the tag assigned to the atom. */
    inline uint_ GetTag() const { return _tag; };
    
    /*! \brief Get number of implicit hydrogens.
     *  \return the number of implicit hydrogens in the atom. */
    inline uint_ GetImplicitCount() const { return _implicitH; }
    
    /*! \brief Add an implicit hydrogen.
     *  \return the new number of implicit hydrogens in the atom. */
    inline uint_ AddImplicitHydrogen() { return ++_implicitH; }
    
    /*! \brief Remove an implicit hydrogen.
     *  \details If GetImplicitCount() == 0, no hydrogen is removed.
     *  \return the new number of implicit hydrogens in the atom. */
    inline uint_ RemoveImplicitHydrogen() {
      return _implicitH ? --_implicitH : _implicitH;
    }
    
    /*! \brief Molecule this atom is associated with.
     *  \details The returned shared_ptr is empty of the atom is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this atom. */
    inline Molecule GetMolecule() const { return _mol.lock(); }
    
    /*! \brief Atom name.
     *  \return name of the atom. */
    inline string_ GetName() const { return _name; }
    
    /*! \brief Atom x position.
     *  \return the x coordinate of this atom. */
    inline float_ GetX() const { return _pos[0]; }
    
    /*! \brief Atom y position.
     *  \return the y coordinate of this atom. */
    inline float_ GetY() const { return _pos[1]; }
    
    /*! \brief Atom z position.
     *  \return the z coordinate of this atom. */
    inline float_ GetZ() const { return _pos[2]; }
    
    /*! \brief Vector of the atom's position.
     *  \return the atoms position. */
    inline const Eigen::Vector3d& GetVector() const { return _pos; }
    
    /*! \brief String representation of the atom.
     *  \details The returned string is of the form: Atom(NAME, SYMBOL).
     *  \return a string representation of the atom. */
    string_ ToString();
    
    /*! \brief Set the element of this atom.
     *  \param e the element to set to. */
    void SetElement(Element e);
    
    /*! \brief Set the element of this atom.
     *  \param e the name or atomic symbol of the element to set. */
    inline void SetElement(string_ e) {
      SetElement(GetPeriodicTable()->GetElement(e));
    }
    
    /*! \brief Set the element of this atom.
     *  \param e the atomic number of the element to set. */
    inline void SetElement(uint_ e){
      SetElement(GetPeriodicTable()->GetElement(e));
    }
    
    /*! \brief Set the formal charge of this atom.
     *  \param q the formal charge value to set. */
    inline void SetFormalCharge(int_ q) { _fc = q; }
    
    /*! \brief Set the partial charge of this atom.
     *  \param q the partial charge value to set. */
    inline void SetPartialCharge(float_ q) { _partial = q; }
    
    /*! \brief Set the number of implicit hydrogens.
     *  \param h the number of implicit hydrogens to set. */
    inline void SetImplicitCount(uint_ h) { _implicitH = h; }
    
    /*! \brief Set the tag of this atom.
     *  \details The tag of an atom should not be considered stable. Use with
     *  caution.
     *  \param i the tag to set. */
    inline void SetTag(uint_ i) { _tag = i; }
    
    /*! \brief Set the atom name.
     *  \param n name to set. */
    inline void SetName(string_ n) { _name = n; }
    
    /*! \brief Set the x position.
     *  \param x position to set. */
    inline void SetX(float_ x) { _pos(0) = x; }
    
    /*! \brief Set the y position.
     *  \param y position to set. */
    inline void SetY(float_ y) { _pos(1) = y; }
    
    /*! \brief Set the z position.
     *  \param z position to set. */
    inline void SetZ(float_ z) { _pos(2) = z; }
    
    /*! \brief Set the x, y and z positions.
     *  \param x,y,z position to set. */
    inline void SetPosition(float_ x, float_ y, float_ z) { _pos << x, y, z; }
    
    /*! \brief Set the stereochemistry of an atomic center.
     *  \param s the stereochemistry to set. */
    inline void SetStereochemistry(Stereo s) { _stereo = s; }
    
    /*! \brief Set the aromaticity of an atom.
     *  \param a if the atom is aromatic or not. */
    inline void SetAromaticity(bool a) { _aromatic = a; }
    
    /*! \brief Get the stereochemistry of the atom.
     *  \return the stereochemistry of the atom. */
    inline Stereo GetStereochemistry() const { return _stereo; }
    
    /*! \brief Get the aromaticity of an atom.
     *  \return if the atom is aromatic or not. */
    inline bool GetAromaticity() const { return _aromatic; }
    
  private:
    /*! \brief Add a bond to this atom.
     *  \details Assumes that the bond is not already added to the atom.
     *  \param b the bond to add.
     *  \return if the bond was added or not. */
    inline void AddBond(Bond b) { _bnds.emplace_back(b); }
    
    /*! \brief Add an angle to this atom.
     *  \details Assumes that the angle is not already added to the atom.
     *  \param a the angle to add.
     *  \return if the angle was added or not. */
    inline void AddAngle(Angle a) { _angs.emplace_back(a); }
    
    /*! \brief Add a dihedral to this atom.
     *  \details Assumes that the dihedral is not already added to the atom.
     *  \param d the dihedral to add.
     *  \return if the dihedral was added or not. */
    inline void AddDihedral(Dihedral d) { _dhds.emplace_back(d); }
    
    /*! \brief Remove a bond from this atom.
     *  \details Assumes that the bond is added to the atom.
     *  \param b the bond to remove.
     *  \return if the bond was removed or not. */
    inline void RemoveBond(Bond b) {
      _bnds.erase(utils::WeakContainsShared(_bnds.begin(), _bnds.end(), b));
    }
    
    /*! \brief Remove an angle from this atom.
     *  \details Assumes that the angle is added to the atom.
     *  \param a the angle to remove.
     *  \return if the angle was removed or not. */
    inline void RemoveAngle(Angle a) {
      _angs.erase(utils::WeakContainsShared(_angs.begin(), _angs.end(), a));
    }
    
    /*! \brief Remove a dihedral from this atom.
     *  \details Assumes that the dihedral is added to the atom.
     *  \param d the dihedral to remove.
     *  \return if the dihedral was removed or not. */
    inline void RemoveDihedral(Dihedral d) {
      _dhds.erase(utils::WeakContainsShared(_dhds.begin(), _dhds.end(), d));
    }
    
    /*! \brief Clear all information.
     *  \details Erases all information stored on the atom. */
    void Clear();
    
  public:
    /*! \brief Get iterator access to the atom's bonds.
     *  \details Intended primarily for internal use as the iterators are to
     *  weak_ptrs.
     *  \returns a pair of iterators for the beginning and end of the bonds. */
    std::pair<AtomBondIter, AtomBondIter> GetBondIters() const {
      return std::make_pair(_bnds.begin(), _bnds.end());
    }
    
    /*! \brief Get iterator access to the atom's angles.
     *  \details Intended primarily for internal use as the iterators are to
     *  weak_ptrs.
     *  \returns a pair of iterators for the beginning and end of the angles. */
    std::pair<AtomAngleIter, AtomAngleIter> GetAngleIters() const {
      return std::make_pair(_angs.begin(), _angs.end());
    }
    
    /*! \brief Get iterator access to the atom's dihedrals.
     *  \details Intended primarily for internal use as the iterators are to
     *  weak_ptrs.
     *  \returns a pair of iterators for the beginning and end of the
     *  dihedrals. */
    std::pair<AtomDihedralIter, AtomDihedralIter> GetDihedralIters() const {
      return std::make_pair(_dhds.begin(), _dhds.end());
    }
    
    /*! \brief Number of valid bonds this atom is part of.
     *  \returns the number of valid assigned bonds. */
    size_ NumBonds() const { return _bnds.size(); }
    
    /*! \brief Number of valid angles this atom is a part of.
     *  \returns the number of valid assigned angles. */
    size_ NumAngles() const { return _angs.size(); }
    
    /*! \brief Number of valid dihedrals this atom is a part of.
     *  \returns the number of valid assigned dihedrals. */
    size_ NumDihedrals() const { return _dhds.size(); }
    
    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the atom in the container of atoms
     *  of the molecule it is a part of. If the molecule is dead, the index
     *  returned is the tag of the atom.
     *  \return the index of the atom. */
    size_ GetIndex() const;
    
    /*! \brief Get the FF type of the atom.
     *  \return the force field type of the atom. */
    FFAtom GetType() const { return _type; }
    
    /*! \brief Set the FF type of the atom.
     *  \param type the type to set. */
    void SetType(FFAtom type) { _type = type; }
    
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
    Eigen::Vector3d _pos;
    //! Partial atomic charge.
    float_ _partial;
    //! Stereochemistry
    Stereo _stereo;
    //! Aromaticity
    bool _aromatic;
    
    //! MM type for atom
    FFAtom _type;
    
    //! Bonds the atom is part of
    AtomBonds _bnds;
    //! Angles the atom is part of
    AtomAngles _angs;
    //! Dihedrals the atom is part of
    AtomDihedrals _dhds;
  };
  
  /*! \brief Print an Atom to an output stream.
   *  \details The printed string is of the form: Atom(UID, SYMBOL).
   *  \param os the output stream to print to.
   *  \param atom the Atom to print.
   *  \return the output stream after printing. */
  inline std::ostream& operator<<(std::ostream& os, const IXAtom& atom) {
    return (os << "Atom(" << atom.GetIndex() << ")");
  }
  inline std::ostream& operator<<(std::ostream& os, const Atom& atom) {
    return atom ? os << *atom : os;
  }
  
  //! \brief Type for the stereochemistry enum of an atom.
  using AtomStereo = indigox::IXAtom::Stereo;
}

#endif /* INDIGOX_CLASSES_ATOM_HPP */
