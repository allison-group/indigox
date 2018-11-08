/*! \file atom.hpp */

#ifndef INDIGOX_CLASSES_ATOM_HPP
#define INDIGOX_CLASSES_ATOM_HPP

#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "forcefield.hpp"
#include "periodictable.hpp"
#include "../utils/common.hpp"
#include "../utils/counter.hpp"
#include "../utils/fwd_declares.hpp"

namespace indigox {
  class Atom
  : public utils::IXCountableObject<Atom>,
  public std::enable_shared_from_this<Atom> {
    //! \brief Friendship allows IXMolecule to create new atoms.
    friend class indigox::Molecule;
    //! \brief Friendship allows IXAtom to be tested.
    friend struct indigox::test::TestAtom;
    //! \brief Friendship allows IXAtom to be serialised.
    friend class cereal::access;
    
  private:
    //! \brief Container for storing Bond references.
    using AtomBonds = std::vector<wBond>;
    //! \brief Container for storing IXAngle references.
    using AtomAngles = std::vector<wAngle>;
    //! \brief Container for storing IXDihedral references.
    using AtomDihedrals = std::vector<wDihedral>;
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
                                   cereal::construct<Atom>& construct,
                                   const uint32_t version);
    
  public:
    Atom() = delete;  // default constructor for serialise access
    Atom(const Atom&) = delete;
    Atom& operator=(const Atom&) = delete;
    /*! \brief Normal constructor.
     *  \details Links the constructed atom to the given Molecule.
     *  \param m the molecule to assign this atom to. */
    Atom(Molecule& m);
    
    //! \brief Destructor
    ~Atom() { };
    
    /*! \brief Element of the atom.
     *  \return the element of this atom. */
    inline const Element& GetElement() const { return _elem; }
    
    /*! \brief Formal charge on the atom.
     *  \return the formal charge on the atom. */
    inline int32_t GetFormalCharge() const { return _fc; }
    
    /*! \brief Partial atomic charge on the atom.
     *  \return the partial atomic charge. */
    inline double GetPartialCharge() const { return _partial; }
    
    /*! \brief Tag of the atom.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the atom, use IXAtom::GetUniqueID.
     *  \return the tag assigned to the atom. */
    inline uint32_t GetTag() const { return _tag; };
    
    /*! \brief Get number of implicit hydrogens.
     *  \return the number of implicit hydrogens in the atom. */
    inline uint32_t GetImplicitCount() const { return _implicitH; }
    
    /*! \brief Add an implicit hydrogen.
     *  \return the new number of implicit hydrogens in the atom. */
    inline uint32_t AddImplicitHydrogen() { return ++_implicitH; }
    
    /*! \brief Remove an implicit hydrogen.
     *  \details If GetImplicitCount() == 0, no hydrogen is removed.
     *  \return the new number of implicit hydrogens in the atom. */
    inline uint32_t RemoveImplicitHydrogen() {
      return _implicitH ? --_implicitH : _implicitH;
    }
    
    /*! \brief Molecule this atom is associated with.
     *  \details The returned shared_ptr is empty of the atom is not assigned
     *  to a valid molecule.
     *  \return the molecule associated with this atom. */
    inline Molecule& GetMolecule() const { return *_mol.lock(); }
    
    /*! \brief Atom name.
     *  \return name of the atom. */
    inline std::string GetName() const { return _name; }
    
    /*! \brief Atom x position.
     *  \return the x coordinate of this atom. */
    inline double GetX() const { return _pos[0]; }
    
    /*! \brief Atom y position.
     *  \return the y coordinate of this atom. */
    inline double GetY() const { return _pos[1]; }
    
    /*! \brief Atom z position.
     *  \return the z coordinate of this atom. */
    inline double GetZ() const { return _pos[2]; }
    
    /*! \brief Vector of the atom's position.
     *  \return the atoms position. */
    inline const Eigen::Vector3d& GetVector() const { return _pos; }
    
    /*! \brief String representation of the atom.
     *  \details The returned string is of the form: Atom(NAME, SYMBOL).
     *  \return a string representation of the atom. */
    std::string ToString();
    
    /*! \brief Set the element of this atom.
     *  \param e the element to set to. */
    void SetElement(const Element& e);
    
    /*! \brief Set the element of this atom.
     *  \param e the name or atomic symbol of the element to set. */
    inline void SetElement(std::string e) {
      SetElement(GetPeriodicTable().GetElement(e));
    }
    
    /*! \brief Set the element of this atom.
     *  \param e the atomic number of the element to set. */
    inline void SetElement(uint32_t e){
      SetElement(GetPeriodicTable().GetElement(e));
    }
    
    /*! \brief Set the formal charge of this atom.
     *  \param q the formal charge value to set. */
    inline void SetFormalCharge(int32_t q) { _fc = q; }
    
    /*! \brief Set the partial charge of this atom.
     *  \param q the partial charge value to set. */
    inline void SetPartialCharge(double q) { _partial = q; }
    
    /*! \brief Set the number of implicit hydrogens.
     *  \param h the number of implicit hydrogens to set. */
    inline void SetImplicitCount(uint32_t h) { _implicitH = h; }
    
    /*! \brief Set the tag of this atom.
     *  \details The tag of an atom should not be considered stable. Use with
     *  caution.
     *  \param i the tag to set. */
    inline void SetTag(uint32_t i) { _tag = i; }
    
    /*! \brief Set the atom name.
     *  \param n name to set. */
    inline void SetName(std::string n) { _name = n; }
    
    /*! \brief Set the x position.
     *  \param x position to set. */
    inline void SetX(double x) { _pos(0) = x; }
    
    /*! \brief Set the y position.
     *  \param y position to set. */
    inline void SetY(double y) { _pos(1) = y; }
    
    /*! \brief Set the z position.
     *  \param z position to set. */
    inline void SetZ(double z) { _pos(2) = z; }
    
    /*! \brief Set the x, y and z positions.
     *  \param x,y,z position to set. */
    inline void SetPosition(double x, double y, double z) { _pos << x, y, z; }
    
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
    void AddBond(Bond& b);
    
    /*! \brief Add an angle to this atom.
     *  \details Assumes that the angle is not already added to the atom.
     *  \param a the angle to add.
     *  \return if the angle was added or not. */
    void AddAngle(Angle& a);
    
    /*! \brief Add a dihedral to this atom.
     *  \details Assumes that the dihedral is not already added to the atom.
     *  \param d the dihedral to add.
     *  \return if the dihedral was added or not. */
    void AddDihedral(Dihedral& d);
    
    /*! \brief Remove a bond from this atom.
     *  \details Assumes that the bond is added to the atom.
     *  \param b the bond to remove.
     *  \return if the bond was removed or not. */
    void RemoveBond(Bond& b);
    
    /*! \brief Remove an angle from this atom.
     *  \details Assumes that the angle is added to the atom.
     *  \param a the angle to remove.
     *  \return if the angle was removed or not. */
    void RemoveAngle(Angle& a);
    
    /*! \brief Remove a dihedral from this atom.
     *  \details Assumes that the dihedral is added to the atom.
     *  \param d the dihedral to remove.
     *  \return if the dihedral was removed or not. */
    void RemoveDihedral(Dihedral& d);
    
    void Clear();
    
  public:
    /*! \brief Get iterator access to the atom's bonds.
     *  \details Intended primarily for internal use as the iterators are to
     *  weak_ptrs.
     *  \returns a pair of iterators for the beginning and end of the bonds. */
    std::pair<AtomBondIter, AtomBondIter> GetBondIters() const {
      return std::make_pair(_bnds.begin(), _bnds.end());
    }
    
    const AtomBonds& GetBonds() const { return _bnds; }
    
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
    size_t NumBonds() const { return _bnds.size(); }
    
    /*! \brief Number of valid angles this atom is a part of.
     *  \returns the number of valid assigned angles. */
    size_t NumAngles() const { return _angs.size(); }
    
    /*! \brief Number of valid dihedrals this atom is a part of.
     *  \returns the number of valid assigned dihedrals. */
    size_t NumDihedrals() const { return _dhds.size(); }
    
    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the atom in the container of atoms
     *  of the molecule it is a part of. If the molecule is dead, the index
     *  returned is the tag of the atom.
     *  \return the index of the atom. */
    size_t GetIndex() const;
    
    /*! \brief Get the FF type of the atom.
     *  \return the force field type of the atom. */
    const FFAtom& GetType() const { return _type; }
    
    bool HasType() const { return bool(_type); }
    
    /*! \brief Set the FF type of the atom.
     *  \param type the type to set. */
    void SetType(const FFAtom& type);
    
  private:
    //! The molecule this atom is assigned to.
    wMolecule _mol;
    //! The atoms element.
    Element _elem;
    //! Formal charge.
    int32_t _fc;
    //! Tag (unstable).
    uint32_t _tag;
    //! Number of implicit hydrogens.
    uint32_t _implicitH;
    //! Atoms name.
    std::string _name;
    //! Position vector.
    Eigen::Vector3d _pos;
    //! Partial atomic charge.
    double _partial;
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
  std::ostream& operator<<(std::ostream& os, const Atom& atom);
  
  //! \brief Type for the stereochemistry enum of an atom.
  using AtomStereo = indigox::Atom::Stereo;
}

#endif /* INDIGOX_CLASSES_ATOM_HPP */
