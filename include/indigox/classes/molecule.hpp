/*! \file molecule.hpp */
#include "../utils/atomic_coordinates.hpp"
#include "../utils/fwd_declares.hpp"
#include "bond.hpp"

#include <bitset>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>
#include <indigo-bondorder/indigo-bondorder.hpp>

#ifndef INDIGOX_CLASSES_MOLECULE_HPP
#define INDIGOX_CLASSES_MOLECULE_HPP

namespace indigox {

  class Molecule {
    //! \brief Friendship allows serialisation
    friend class cereal::access;

  public:
    /*! \brief Container for storing IXAtom instances.
     *  \details A molecule takes ownership of all atoms it contains. */
    using MoleculeAtoms = std::vector<Atom>;
    /*! \brief Container for storing IXBond instances.
     *  \details A molecule takes ownership of all bonds it contains. */
    using MoleculeBonds = std::vector<Bond>;
    /*! \brief Container for storing IXAngle instances.
     *  \details A molecule takes ownership of all angles it contains. */
    using MoleculeAngles = std::vector<Angle>;
    /*! \brief Container for storing IXDihedral instances.
     *  \details A molecule takes ownership of all dihedrals it contains. */
    using MoleculeDihedrals = std::vector<Dihedral>;

    using MoleculeResidues = std::vector<Residue>;

  private:
    template <typename Archive>
    void serialise(Archive &archive, const uint32_t version);

  public:
    /*! \brief Default constructor
     *  \details Is private to enforce that IXMolecules should only be used
     *  via the Molecule shared_ptr. */
    Molecule(const std::string& name);

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Molecule);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Molecule, mol);

  public:
    /*! \brief Check if the atom is owned by this molecule.
     *  \param atom the atom to check for.
     *  \return if the atom is owned by this molecule. */
    bool HasAtom(const Atom &atom) const;

    /*! \brief Check if the bond is owned by this molecule.
     *  \param bond the bond to check for.
     *  \return if the bond is owned by this molecule. */
    bool HasBond(const Bond &bond) const;

    /*! \brief Check if a bond between two atoms exists in this molecule.
     *  \param a,b the atoms the bond should be between.
     *  \return if there is a bond between the two atoms. */
    bool HasBond(const Atom &a, const Atom &b) const;

    /*! \brief Check if the angle is owned by this molecule.
     *  \param angle the angle to check for.
     *  \return if the angle is owned by this molecule. */
    bool HasAngle(const Angle &angle) const;

    /*! \brief Check if an angle between three atoms exists in this molecule.
     *  \details If the Emergent::ANGLE_PERCEPTION property is set, angles are
     *  re-percieved prior to check for the existance of the angle. \p b is the
     *  central atom of the angle.
     *  \param a,b,c the atoms the angle should be between.
     *  \return if there is an angle between the three atoms. */
    bool HasAngle(const Atom &a, const Atom &b, const Atom &c);

    /*! \brief Check if the dihedral is owned by this molecule.
     *  \param dihedral the dihedral to check for.
     *  \return if the dihedral is owned by this molecule. */
    bool HasDihedral(const Dihedral &dihedral) const;

    /*! \brief Check if a dihedral between four atoms exists in this molecule.
     *  \details If the Emergent::DIHEDRAL_PERCEPTION property is set, dihedrals
     *  are re-perceived prior to checking for the existance of the dihedral.
     *  \param a,b,c,d the atoms the dihedral should be between.
     *  \return if there is a dihedral between the four atoms. */
    bool HasDihedral(const Atom &a, const Atom &b, const Atom &c,
                     const Atom &d);

    /*! \brief Get the number of atoms in the molecule.
     *  \return the number of atoms in the molecule. */
    int64_t NumAtoms() const;

    /*! \brief Get the number of bonds in the molecule.
     *  \return the number of bonds in the molecule. */
    int64_t NumBonds() const;

    /*! \brief Get the number of angles in the molecule.
     *  \details If the Emergent::ANGLE_PERCEPTION property is set, angles are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of angles in the molecule. */
    int64_t NumAngles();

    /*! \brief Get the number of dihedrals in the molecule.
     *  \details If the CONNECTIVITY property has been modified, dihedrals are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of dihedrals in the molecule. */
    int64_t NumDihedrals();

  public:
    /*! \brief Get the atom at position \p pos.
     *  \details Returns that atom at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of atoms may
     *  change during normal operations.
     *  \param pos the position of the atom to get.
     *  \return the atom at pos or an empty shared_ptr. */
    Atom GetAtom(uint32_t pos) const;

    /*! \brief Get the atom with the given id.
     *  \details Returns the atom with the given unique id. If no such atom is
     *  found returns an empty shared_ptr.
     *  \param id the unique id of the atom to retrieve.
     *  \return the atom with the given id or an empty shared_ptr. */
    Atom GetAtomID(int64_t id) const;

    /*! \brief Get the first atom with the given tag.
     *  \details Returns the first atom with a tag matching that given. If no
     *  such atom is found, returns an empty shared_ptr.
     *  \param tag the atom tag to search for.
     *  \return the first atom with the tag or an empty shared_ptr. */
    Atom GetAtomTag(int64_t tag) const;

    /*! \brief Get the bond at position \p pos.
     *  \details Returns that bond at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of bonds may
     *  change during normal operations.
     *  \param pos the position of the bond to get.
     *  \return the bond at pos or an empty shared_ptr. */
    Bond GetBond(uint32_t pos) const;

    /*! \brief Get the bond between two atoms, if it exists.
     *  \details Searches for a bond between \p a and \p b and returns it. If no
     *  bond is found, returns an empty shared_ptr.
     *  \param a,b the atoms to get the bond between.
     *  \return the bond between \p a and \p b or an empty shared_ptr. */
    Bond GetBond(const Atom &a, const Atom &b) const;

    /*! \brief Get the bond with the given id.
     *  \details Returns the bond with the given unique id. If no such bond is
     *  found returns an empty shared_ptr.
     *  \param id the unique id of the bond to retrieve.
     *  \return the bond with the given id or an empty shared_ptr. */
    Bond GetBondID(int64_t id) const;

    /*! \brief Get the first bond with the given tag.
     *  \details Returns the first bond with a tag matching that given. If no
     *  such bond is found, returns an empty shared_ptr.
     *  \param tag the bond tag to search for.
     *  \return the first bond with the tag or an empty shared_ptr. */
    Bond GetBondTag(int64_t tag) const;

    /*! \brief Get the angle at position \p pos.
     *  \details Returns the angle at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of angles
     * may change during normal operations. Angles are not re-percieved prior to
     *  retrieving the angle.
     *  \param pos the position of the angle to get.
     *  \return the angle at pos or an empty shared_ptr. */
    Angle GetAngle(uint32_t pos);

    /*! \brief Get the angle between three atoms, if it exists.
     *  \details Searches for an angle between the three atoms and returns it.
     *  If no angle is found, the returned shared_ptr is empty. Angles are
     *  re-perceived prior to searching.
     *  \param a,b,c the atoms to get the angle between. \p b is the central
     *  atom.
     *  \return the angle between \p a and \p c centred on \p b. */
    Angle GetAngle(const Atom &a, const Atom &b, const Atom &c);

    /*! \brief Get the angle with the given id.
     *  \details Returns the angle with the given unique id. If no such angle is
     *  found returns an empty shared_ptr. As the unique id of any created
     * angles will be unknown, angles are not re-percieved prior to searching.
     *  \param id the unique id of the angle to retrieve.
     *  \return the angle with the given id or an empty shared_ptr. */
    Angle GetAngleID(int64_t id) const;

    /*! \brief Get the first angle with the given tag.
     *  \details Returns the first angle with a tag matching that given. If no
     *  such angle is found, returns an empty shared_ptr. As tags are not set on
     *  creation of angles, angles are not re-percieved prior to searching.
     *  \param tag the angle tag to search for.
     *  \return the first angle with the tag or an empty shared_ptr. */
    Angle GetAngleTag(int64_t tag) const;

    /*! \brief Get the angle at position \p pos.
     *  \details Returns the dihedral at \p pos after a range check. If \p pos
     *  is not a valid index, returns an empty shared_ptr. Positioning of
     *  dihedrals may change during normal operations. Dihedrals are not
     *  re-perceived prior to retrieving the dihedral.
     *  \param pos the position of the dihedral to get.
     *  \return the dihedral at pos or an empty shared_ptr. */
    Dihedral GetDihedral(uint32_t pos);

    /*! \brief Get the dihedral between four atoms, it if exists.
     *  \details Searches for a dihedral between the four atoms and returns it.
     *  If no dihedral is found, the returned shared_ptr is empty. Dihedrals are
     *  re-perceived prior to searching.
     *  \param a,b,c,d the atoms to get the dihedral between.
     *  \return the dihedral between \p a, \p b, \p c and \p d. */
    Dihedral GetDihedral(const Atom &a, const Atom &b, const Atom &c,
                         const Atom &d);

    /*! \brief Get the dihedral with the given id.
     *  \details Returns the dihedral with the given unique id. If no such
     *  dihedral is found returns an empty shared_ptr. As the unique id of any
     *  created dihedrals will be unknown, dihedrals are not re-perceived prior
     *  to searching.
     *  \param id the unique id of the dihedral to retrieve.
     *  \return the dihedral with the given id or an empty shared_ptr. */
    Dihedral GetDihedralID(int64_t id) const;

    /*! \brief Get the first dihedral with the given tag.
     *  \details Returns the first dihedral with a tag matching that given. If
     *  no such dihedral is found, returns an empty shared_ptr. As tags are not
     *  set on creation of dihedrals, dihedrals are not re-perceived prior to
     *  searching.
     *  \param tag the dihedral tag to search for.
     *  \return the first dihedral with the tag or an empty shraed_ptr. */
    Dihedral GetDihedralTag(int64_t tag) const;

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
    const graph::MolecularGraph &GetGraph() const;

    const graph::CondensedMolecularGraph &GetCondensedGraph() const;

    /*! \brief Get the name of the molecule.
     *  \return the name of the molecule. */
    const std::string &GetName() const;

    /*! \brief Get the molecular charge of the molecule.
     *  \return the molecular charge of the molecule. */
    int32_t GetMolecularCharge() const;

    /*! \brief Set the name of the molecule.
     *  \param name the new to set. */
    void SetName(std::string name);

    /*! \brief Set the molecular charge of the molecule.
     *  \details Sets the IXMolecule::ELECTRON_COUNT property as modified.
     *  \param q the new charge to set. */
    void SetMolecularCharge(int32_t q);

    /*! \brief Create a new atom owned by the molecule.
     *  \return the new atom. */
    Atom NewAtom();

    /*! \brief Create a new atom of the given element for the molecule.
     *  \param element the element of the new atom.
     *  \return the new atom. */
    Atom NewAtom(const Element &element);

    /*! \brief Create a new named atom of the given element for the molecule.
     *  \param name the name of the new atom.
     *  \param element the element of the new atom.
     *  \return the new atom. */
    Atom NewAtom(const Element &element, double x, double y, double z);

    /*! \brief Create a bond between two atoms.
     *  \details To create a bond, both atoms need to be owned by the molecule
     *  and there cannot be an existing bond between the atoms. If either
     *  condition is not met, the returned shared_ptr is empty.
     *  \param a,b the atoms to create a bond between.
     *  \return the new bond. */
    Bond NewBond(const Atom &a, const Atom &b);

  private:
    /*! \brief Create an angle between three atoms.
     *  \param a,b,c the atoms to create the angle between.
     *  \return the new angle. */
    Angle NewAngle(const Atom &a, const Atom &b, const Atom &c);

    Dihedral NewDihedral(const Atom &a, const Atom &b, const Atom &c,
                         const Atom &d, bool manual);

  public:
    /*! \brief Create a dihedral between four atoms.
     *  \param a,b,c,d the atoms to create the dihedral between.
     *  \return the new dihedral. */
    Dihedral NewDihedral(const Atom &a, const Atom &b, const Atom &c,
                         const Atom &d);

    /*! \brief Remove an atom from the molecule.
     *  \details In addition to removing the atom, any bonds to the atom are
     *  also removed. A removal will not occur if the atom is not owned by the
     *  molecule it is being removed from.
     *  \param atom the atom to remove.
     *  \return if a removal occured. */
    bool RemoveAtom(const Atom &atom);

    /*! \brief Remove a bond from the molecule.
     *  \details A removal will not occur if the bond is not owned by the
     *  molecule it is being removed from. Reference to the bond is removed from
     *  both of the atoms of the bond.
     *  \param bond the bond to remove.
     *  \return if a removal occured. */
    bool RemoveBond(const Bond &bond);

    /*! \brief Remove a bond between two atoms.
     *  \param a,b the atoms to remove a bond between.
     *  \return if removal occured. */
    bool RemoveBond(const Atom &a, const Atom &b);

  public:
    /*! \brief Determine angles in the molecule.
     *  \details An angle is defined for each group of three connected atoms
     *  within a molecule. Subsequent calls to this method will only generate
     *  angles which have not been previously generated. Additionally, angles
     *  can only be removed by removing an atom or a bond.
     *  \return the number of angles added. */
    int64_t PerceiveAngles();

    /*! \brief Determine dihedrals in the molecule.
     *  \details A dihedral is defined for each group of four connected atoms
     *  within a molecule. Subsequent calls to this method will only generate
     *  dihedrals which have not been previously generated. Additionally,
     *  dihedrals can only be removed by removing an atom or bond.
     *  \return the number of dihedrals added. */
    int64_t PerceiveDihedrals();

    /*! \brief Determine bond order and formal charges in the molecule.
     *  \details Assign electrons to the molecule and determine bond order and
     *  formal charges. Uses the algorithm described at https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0340-0
     *  \param algorithmOption Integer option representing the algorithm to use.
     *  0 = Local Optimisation, 1 = A*, 2 = FPT (fixed parameter tractable)
     *  \return the number of resonance structures found. */
    int64_t PerceiveElectrons(int32_t algorithmOption, bool silent);

    int32_t PerceiveResidues();
    //    size_t AssignElectrons();
    //    bool ApplyElectronAssignment(size_t);
    //    FCSCORE GetMinimumElectronAssignmentScore();

    /*! \brief Reserve storage space for atoms.
     *  \details Reserves storage space for a minimum of \p num IXAtom
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more atoms are added.
     *  \param num the number of IXAtoms to reserve space for. */
    void ReserveAtoms(int64_t num);

    /*! \brief Reserve storage space for bonds.
     *  \details Reserves storage space for a minimum of \p num IXBond
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more bonds are added.
     *  \param num the number of IXBonds to reserve space for. */
    void ReserveBonds(int64_t num);

    const MoleculeAtoms &GetAtoms() const;

    const MoleculeBonds &GetBonds() const;

    const MoleculeAngles &GetAngles(); // percevie first

    const MoleculeDihedrals &GetDihedrals(); // perceive first

    Residue GetResidueID(int32_t id); // perceive first

    const MoleculeResidues &GetResidues(); // perceive first

    const Forcefield &GetForcefield() const;

    void SetForcefield(const Forcefield &ff);

    void ResetForcefield(const Forcefield &ff);

    bool HasForcefield() const;

    void OptimiseChargeGroups();

    void ReorderAtoms(MoleculeAtoms &new_order);

    void ModificationMade();

    AtomicCoordinates &GetAtomicCoordinates();
    
    void UniquifyAtomNames();
    
    void GiveAromaticBondsImpropers();

  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;

    static std::map<indigo_bondorder::BondOrder, indigox::Bond::Order> GetBondorderEnumMap();
    static std::map<indigox::Bond::Order, std::string> GetBondorderNameMap();

    /*!
     * Translate CherryPicker settings into indigo-bondorder settings for the electron calculation
     */
    static void setElectronSettings(int32_t algorithmOption);

    uint chooseResonanceStructure(indigo_bondorder::Molecule_p &mol,
                                  const std::map<indigo_bondorder::BondOrder, indigox::Bond::Order> &enum_map,
                                  const std::map<indigox::Bond::Order, std::string> &name_map,
                                  indigo_bondorder::Uint num_structures);

    static void displayResonanceStructures(const indigo_bondorder::Molecule_p &mol, indigo_bondorder::Uint num_structures,
                                              std::map<indigo_bondorder::BondOrder, indigox::Bond::Order> enum_map,
                                              std::map<indigox::Bond::Order, std::string> string_map);

    static int32_t getChoiceOfStructure();

    static std::string trimOrFill(std::string str, int length);

    static std::string getNameAndIndex(const std::shared_ptr<indigo_bondorder::Atom> &atom);
  };

  void SaveMolecule(const Molecule &mol, std::string path);
  Molecule LoadMolecule(std::string path);

} // namespace indigox

#endif /* INDIGOX_CLASSES_MOLECULE_HPP */
