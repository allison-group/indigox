/*! \file molecule.hpp */
#include <bitset>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

#include "angle.hpp"
#include "atom.hpp"
#include "bond.hpp"
#include "dihedral.hpp"
#include "periodictable.hpp"
#include "../graph/molecular.hpp"
#include "../utils/counter.hpp"
#include "../utils/fwd_declares.hpp"
#include "../utils/modifable_object.hpp"

#ifndef INDIGOX_MOLECULE_HPP
#define INDIGOX_MOLECULE_HPP

namespace indigox {
 
  class Molecule :
  public utils::IXCountableObject<Molecule>,
  public utils::ModifiableObject,
  public std::enable_shared_from_this<Molecule> {
    //! \brief Friendship allows generation of molecules.
    friend sMolecule CreateMolecule();
    //! \brief Friendship allows IXMolecule internals to be tested.
    friend struct indigox::test::TestMolecule;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
  public:
    /*! \brief Container for storing IXAtom instances.
     *  \details A molecule takes ownership of all atoms it contains. */
    using MolAtoms = std::vector<sAtom>;
    /*! \brief Container for storing IXBond instances.
     *  \details A molecule takes ownership of all bonds it contains. */
    using MolBonds = std::vector<sBond>;
    /*! \brief Container for storing IXAngle instances.
     *  \details A molecule takes ownership of all angles it contains. */
    using MolAngles = std::vector<sAngle>;
    /*! \brief Container for storing IXDihedral instances.
     *  \details A molecule takes ownership of all dihedrals it contains. */
    using MolDihedrals = std::vector<sDihedral>;
    
  public:   // Public iterator aliases for easier external usage
    //! \brief Iterator over owned IXAtom instances.
    using MolAtomIter = MolAtoms::const_iterator;
    //! \brief Iterator over owned IXBond instances.
    using MolBondIter = MolBonds::const_iterator;
    //! \brief Iterator over owned IXAngle instances.
    using MolAngleIter = MolAngles::const_iterator;
    //! \brief Iterator over owned IXDihedral instances.
    using MolDihedralIter = MolDihedrals::const_iterator;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    void load(Archive& archive, const uint32_t version);
    
  private:
    /*! \brief Default constructor
     *  \details Is private to enforce that IXMolecules should only be used
     *  via the Molecule shared_ptr. */
    Molecule();
    
  public:
    //! \brief Destructor to clear all members, just in case
    ~Molecule();
    
    /*! \brief Get the angle at position \p pos.
     *  \details Returns the angle at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of angles may
     *  change during normal operations. Angles are not re-percieved prior to
     *  retrieving the angle.
     *  \param pos the position of the angle to get.
     *  \return the angle at pos or an empty shared_ptr. */
    inline Angle& GetAngle(size_t pos) const { return *_angs[pos]; }
    
    /*! \brief Get the angle at position \p pos.
     *  \details Returns the dihedral at \p pos after a range check. If \p pos
     *  is not a valid index, returns an empty shared_ptr. Positioning of
     *  dihedrals may change during normal operations. Dihedrals are not
     *  re-perceived prior to retrieving the dihedral.
     *  \param pos the position of the dihedral to get.
     *  \return the dihedral at pos or an empty shared_ptr. */
    inline Dihedral& GetDihedral(size_t pos) const { return *_dhds[pos]; }
    
    /*! \brief Get the angle between three atoms, if it exists.
     *  \details Searches for an angle between the three atoms and returns it.
     *  If no angle is found, the returned shared_ptr is empty. Angles are
     *  re-perceived prior to searching.
     *  \param a,b,c the atoms to get the angle between. \p b is the central
     *  atom.
     *  \return the angle between \p a and \p c centred on \p b. */
    Angle& GetAngle(Atom& a, Atom& b, Atom& c);
    
    /*! \brief Get the dihedral between four atoms, it if exists.
     *  \details Searches for a dihedral between the four atoms and returns it.
     *  If no dihedral is found, the returned shared_ptr is empty. Dihedrals are
     *  re-perceived prior to searching.
     *  \param a,b,c,d the atoms to get the dihedral between.
     *  \return the dihedral between \p a, \p b, \p c and \p d. */
    Dihedral& GetDihedral(Atom& a, Atom& b, Atom& c, Atom& d);
    
    /*! \brief Get the first angle with the given tag.
     *  \details Returns the first angle with a tag matching that given. If no
     *  such angle is found, returns an empty shared_ptr. As tags are not set on
     *  creation of angles, angles are not re-percieved prior to searching.
     *  \param tag the angle tag to search for.
     *  \return the first angle with the tag or an empty shared_ptr. */
    Angle& GetAngleTag(uint32_t tag) const;
    
    /*! \brief Get the first dihedral with the given tag.
     *  \details Returns the first dihedral with a tag matching that given. If
     *  no such dihedral is found, returns an empty shared_ptr. As tags are not
     *  set on creation of dihedrals, dihedrals are not re-perceived prior to
     *  searching.
     *  \param tag the dihedral tag to search for.
     *  \return the first dihedral with the tag or an empty shraed_ptr. */
    Dihedral& GetDihedralTag(uint32_t tag) const;
    
    /*! \brief Get the angle with the given id.
     *  \details Returns the angle with the given unique id. If no such angle is
     *  found returns an empty shared_ptr. As the unique id of any created angles
     *  will be unknown, angles are not re-percieved prior to searching.
     *  \param id the unique id of the angle to retrieve.
     *  \return the angle with the given id or an empty shared_ptr. */
    Angle& GetAngleID(uint32_t id) const;
    
    /*! \brief Get the dihedral with the given id.
     *  \details Returns the dihedral with the given unique id. If no such
     *  dihedral is found returns an empty shared_ptr. As the unique id of any
     *  created dihedrals will be unknown, dihedrals are not re-perceived prior
     *  to searching.
     *  \param id the unique id of the dihedral to retrieve.
     *  \return the dihedral with the given id or an empty shared_ptr. */
    Dihedral& GetDihedralID(uint32_t id) const;
    
    /*! \brief Get the atom at position \p pos.
     *  \details Returns that atom at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of atoms may
     *  change during normal operations.
     *  \param pos the position of the atom to get.
     *  \return the atom at pos or an empty shared_ptr. */
    inline Atom& GetAtom(size_t pos) const { return *_atms[pos]; }
    
    /*! \brief Get the first atom with the given tag.
     *  \details Returns the first atom with a tag matching that given. If no
     *  such atom is found, returns an empty shared_ptr.
     *  \param tag the atom tag to search for.
     *  \return the first atom with the tag or an empty shared_ptr. */
    Atom& GetAtomTag(uint32_t tag) const;
    
    /*! \brief Get the atom with the given id.
     *  \details Returns the atom with the given unique id. If no such atom is
     *  found returns an empty shared_ptr.
     *  \param id the unique id of the atom to retrieve.
     *  \return the atom with the given id or an empty shared_ptr. */
    Atom& GetAtomID(uint32_t id) const;
    
    /*! \brief Get the bond at position \p pos.
     *  \details Returns that bond at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of bonds may
     *  change during normal operations.
     *  \param pos the position of the bond to get.
     *  \return the bond at pos or an empty shared_ptr. */
    inline Bond& GetBond(size_t pos) const { return *_bnds[pos]; }
    
    /*! \brief Get the bond between two atoms, if it exists.
     *  \details Searches for a bond between \p a and \p b and returns it. If no
     *  bond is found, returns an empty shared_ptr.
     *  \param a,b the atoms to get the bond between.
     *  \return the bond between \p a and \p b or an empty shared_ptr. */
    Bond& GetBond(Atom& a, Atom& b) const;
    
    /*! \brief Get the first bond with the given tag.
     *  \details Returns the first bond with a tag matching that given. If no
     *  such bond is found, returns an empty shared_ptr.
     *  \param tag the bond tag to search for.
     *  \return the first bond with the tag or an empty shared_ptr. */
    Bond& GetBondTag(uint32_t tag) const;
    
    /*! \brief Get the bond with the given id.
     *  \details Returns the bond with the given unique id. If no such bond is
     *  found returns an empty shared_ptr.
     *  \param id the unique id of the bond to retrieve.
     *  \return the bond with the given id or an empty shared_ptr. */
    Bond& GetBondID(uint32_t id) const;
    
    /*! \brief Get the molecular formula of the molecule.
     *  \details Determines the molecular formula of the molecule. The elements
     *  of the formula are arranged in the order: carbon, hydrogen and then
     *  alphabetically by atomic symbol. The resultant formula is cached and
     *  recalculated only when the IXMolecule::ATOM_ELEMENTS property has been
     *  modified.
     *  \return the molecular formula of the molecule. */
    std::string GetFormula();
    
    /*! \brief Get the molecular graph for this molecule.
     *  \return the molecular graph of this molecule. */
    inline const graph::MolecularGraph& GetGraph() const { return *_g; }
    
    /*! \brief Get the name of the molecule.
     *  \return the name of the molecule. */
    inline std::string GetName() const { return _name; }
    
    /*! \brief Get the molecular charge of the molecule.
     *  \return the molecular charge of the molecule. */
    inline int GetMolecularCharge() const { return _q; }
    
    /*! \brief Get the number of atoms in the molecule.
     *  \return the number of atoms in the molecule. */
    inline size_t NumAtoms() const { return _atms.size(); }
    
    /*! \brief Get the number of bonds in the molecule.
     *  \return the number of bonds in the molecule. */
    inline size_t NumBonds() const { return _bnds.size(); }
    
    /*! \brief Get the number of angles in the molecule.
     *  \details If the Emergent::ANGLE_PERCEPTION property is set, angles are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of angles in the molecule. */
    inline size_t NumAngles() { PerceiveAngles(); return _angs.size(); }
    
    /*! \brief Get the number of dihedrals in the molecule.
     *  \details If the CONNECTIVITY property has been modified, dihedrals are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of dihedrals in the molecule. */
    inline size_t NumDihedrals() { PerceiveDihedrals(); return _dhds.size(); }
    
    /*! \brief Set the name of the molecule.
     *  \param name the new to set. */
    inline void SetName(std::string name) { _name = name; }
    
    /*! \brief Set the molecular charge of the molecule.
     *  \details Sets the IXMolecule::ELECTRON_COUNT property as modified.
     *  \param q the new charge to set. */
    inline void SetMolecularCharge(int q) { _q = q; }
    
    /*! \brief Check if the atom is owned by this molecule.
     *  \param atom the atom to check for.
     *  \return if the atom is owned by this molecule. */
    inline bool HasAtom(Atom& atom) const {
      return &atom.GetMolecule() == this; }
    
    /*! \brief Check if the bond is owned by this molecule.
     *  \param bond the bond to check for.
     *  \return if the bond is owned by this molecule. */
    inline bool HasBond(Bond& bond) const {
      return &bond.GetMolecule() == this; }
    
    /*! \brief Check if a bond between two atoms exists in this molecule.
     *  \param a,b the atoms the bond should be between.
     *  \return if there is a bond between the two atoms. */
    bool HasBond(Atom& a, Atom& b) const;
    
    /*! \brief Check if the angle is owned by this molecule.
     *  \param angle the angle to check for.
     *  \return if the angle is owned by this molecule. */
    inline bool HasAngle(Angle& angle) const {
      return &angle.GetMolecule() == this;
    }
    
    /*! \brief Check if an angle between three atoms exists in this molecule.
     *  \details If the Emergent::ANGLE_PERCEPTION property is set, angles are
     *  re-percieved prior to check for the existance of the angle. \p b is the
     *  central atom of the angle.
     *  \param a,b,c the atoms the angle should be between.
     *  \return if there is an angle between the three atoms. */
    bool HasAngle(Atom& a, Atom& b, Atom& c);
    
    /*! \brief Check if the dihedral is owned by this molecule.
     *  \param dihedral the dihedral to check for.
     *  \return if the dihedral is owned by this molecule. */
    inline bool HasDihedral(Dihedral& dihedral) const {
      return &dihedral.GetMolecule() == this;
    }
    
    /*! \brief Check if a dihedral between four atoms exists in this molecule.
     *  \details If the Emergent::DIHEDRAL_PERCEPTION property is set, dihedrals
     *  are re-perceived prior to checking for the existance of the dihedral.
     *  \param a,b,c,d the atoms the dihedral should be between.
     *  \return if there is a dihedral between the four atoms. */
    bool HasDihedral(Atom& a, Atom& b, Atom& c, Atom& d);
    
    /*! \brief Create a new atom owned by the molecule.
     *  \return the new atom. */
    Atom& NewAtom();
    
    /*! \brief Create a new atom of the given element for the molecule.
     *  \param element the element of the new atom.
     *  \return the new atom. */
    Atom& NewAtom(Element& element);
    
    /*! \brief Create a new named atom owned by the molecule.
     *  \param name the name of the new atom.
     *  \return the new atom. */
    Atom& NewAtom(std::string name);
    
    /*! \brief Create a new named atom of the given element for the molecule.
     *  \param name the name of the new atom.
     *  \param element the element of the new atom.
     *  \return the new atom. */
    Atom& NewAtom(std::string name, Element& element);
    
    /*! \brief Create a bond between two atoms.
     *  \details To create a bond, both atoms need to be owned by the molecule
     *  and there cannot be an existing bond between the atoms. If either
     *  condition is not met, the returned shared_ptr is empty.
     *  \param a,b the atoms to create a bond between.
     *  \return the new bond. */
    Bond& NewBond(Atom& a, Atom& b);
    
    /*! \brief Remove an atom from the molecule.
     *  \details In addition to removing the atom, any bonds to the atom are
     *  also removed. A removal will not occur if the atom is not owned by the
     *  molecule it is being removed from.
     *  \param atom the atom to remove.
     *  \return if a removal occured. */
    bool RemoveAtom(Atom& atom);
    
    /*! \brief Remove a bond from the molecule.
     *  \details A removal will not occur if the bond is not owned by the
     *  molecule it is being removed from. Reference to the bond is removed from
     *  both of the atoms of the bond.
     *  \param bond the bond to remove.
     *  \return if a removal occured. */
    bool RemoveBond(Bond& bond);
    
    /*! \brief Remove a bond between two atoms.
     *  \param a,b the atoms to remove a bond between.
     *  \return if removal occured. */
    inline bool RemoveBond(Atom& a, Atom& b) { return RemoveBond(GetBond(a,b)); }
    
  private:
    /*! \brief Create an angle between three atoms.
     *  \param a,b,c the atoms to create the angle between.
     *  \return the new angle. */
    void NewAngle(Atom& a, Atom& b, Atom& c);
    
  public:
    /*! \brief Create a dihedral between four atoms.
     *  \param a,b,c,d the atoms to create the dihedral between.
     *  \return the new dihedral. */
    Dihedral& NewDihedral(Atom& a, Atom& b, Atom& c, Atom& d);
    
  private:
    //! \cond
    // Return the index of found, -1 if not found
    int64_t _FindBond(Atom& a, Atom& b) const;
    int64_t _FindAngle(Atom& a, Atom& b, Atom& c) const;
    int64_t _FindDihedral(Atom& a, Atom& b, Atom& c, Atom& d) const;
    //! \endcond
    
  public:
    /*! \brief Determine angles in the molecule.
     *  \details An angle is defined for each group of three connected atoms
     *  within a molecule. Subsequent calls to this method will only generate
     *  angles which have not been previously generated. Additionally, angles
     *  can only be removed by removing an atom or a bond.
     *  \return the number of angles added. */
    size_t PerceiveAngles();
    
    /*! \brief Determine dihedrals in the molecule.
     *  \details A dihedral is defined for each group of four connected atoms
     *  within a molecule. Subsequent calls to this method will only generate
     *  dihedrals which have not been previously generated. Additionally,
     *  dihedrals can only be removed by removing an atom or bond.
     *  \return the number of dihedrals added. */
    size_t PerceiveDihedrals();
//    size_t AssignElectrons();
//    bool ApplyElectronAssignment(size_t);
//    FCSCORE GetMinimumElectronAssignmentScore();
    
    /*! \brief Reserve storage space for atoms.
     *  \details Reserves storage space for a minimum of \p num IXAtom
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more atoms are added.
     *  \param num the number of IXAtoms to reserve space for. */
    inline void ReserveAtoms(size_t num) {
      if (_atms.size() < num) _atms.reserve(num);
    }
    
    /*! \brief Reserve storage space for bonds.
     *  \details Reserves storage space for a minimum of \p num IXBond
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more bonds are added.
     *  \param num the number of IXBonds to reserve space for. */
    inline void ReserveBonds(size_t num) {
      if (_bnds.size() < num ) _bnds.reserve(num);
    }
    
    /*! \brief Get iterator access to the owned atoms.
     *  \return a pair of iterators indicating the begining and end of the
     *  owned atoms. */
    inline std::pair<MolAtomIter, MolAtomIter> GetAtomIters() const {
      return std::make_pair(_atms.begin(), _atms.end());
    }
    
    inline const MolAtoms& GetAtoms() const { return _atms; }
    
    /*! \brief Get iterator access to the owned bonds.
     *  \return a pair of iterators indication the beginning and end of the
     *  owned bonds. */
    inline std::pair<MolBondIter, MolBondIter> GetBondIters() const {
      return std::make_pair(_bnds.begin(), _bnds.end());
    }
    
    inline const MolBonds& GetBonds() const { return _bnds; }
    
    /*! \brief Get iterator access to the owned angles.
     *  \details Angles are re-percieved prior to returning the iterators.
     *  \return a pair of iterators for the beginning and end of the angles. */
    inline std::pair<MolAngleIter, MolAngleIter> GetAngleIters() {
      PerceiveAngles();
      return std::make_pair(_angs.begin(), _angs.end());
    }
    
    inline const MolAngles& GetAngles() const { return _angs; }
    
    /*! \brief Get iterator access to the owned dihedrals.
     *  \details Dihedrals are re-perceived prior to returning the iterators.
     *  \return a pair of iterators indication the beginning and end of the
     *  owned dihedrals. */
    inline std::pair<MolDihedralIter, MolDihedralIter> GetDihedralIters() {
      PerceiveDihedrals();
      return std::make_pair(_dhds.begin(), _dhds.end());
    }
    
    inline const MolDihedrals& GetDihedrals() const { return _dhds; }
    
  private:
    //! Name of the molecule
    std::string _name;
    //! Molecular charge of the molecule
    int _q;
    //! Atoms owned by molecule
    MolAtoms _atms;
    //! Bonds owned by molecule
    MolBonds _bnds;
    //! Angles owned by molecule
    MolAngles _angs;
    //! Dihedrals owned by molecule
    MolDihedrals _dhds;
    //! Molecular graph of molecule
    graph::sMolecularGraph _g;
    //! Cached molecular formula string
    std::pair<ModifiableObject::State, std::string> _formula_cache;
    //! State when angles were last percieved
    ModifiableObject::State _angle_perceive;
    //! State when dihedrals were last percieved
    ModifiableObject::State _dihedral_perceive;
  };
}

#endif /* INDIGOX_MOLECULE_HPP */
