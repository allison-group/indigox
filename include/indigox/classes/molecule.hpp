/*! \file molecule.hpp */
#include <bitset>
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
#include "../utils/numerics.hpp"

#ifndef INDIGOX_MOLECULE_HPP
#define INDIGOX_MOLECULE_HPP

namespace indigox {

  class IXMolecule;
  namespace test { struct TestMolecule; }

  //! \brief shared_ptr for normal use of the IXMolecule class.
  using Molecule = std::shared_ptr<IXMolecule>;

  /*! \brief weak_ptr for non-ownership reference to the IXMolecule class.
   *  \details Intended for internal use only. */
  using _Molecule = std::weak_ptr<IXMolecule>;
  
  class IXMolecule : public utils::IXCountableObject<IXMolecule>,
  public std::enable_shared_from_this<IXMolecule> {
    //! \brief Friendship allows generation of molecules.
    friend Molecule CreateMolecule();
    //! \brief Friendship allows IXMolecule internals to be tested.
    friend struct indigox::test::TestMolecule;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
  public: // Public so that IXAtom etc can set when they're modified.
    /*! \brief Enum for the different types of properties a molecule has.
     *  \details A #Property is a minor aspect of a molecule which, if modified,
     *  will potentially require recalculation of #Emergent properties. Each
     *  #Property causes different sets of #Emergent properties to need to be
     *  recalculated. */
    enum class Property : size_t {
      ATOM_ELEMENTS,    //!< Used when elements are changed.
      CONNECTIVITY,     //!< Used when bonds are added/removed.
      ELECTRON_COUNT,   //!< Used when the number of electrons has changed.
      NUM_PROPERTIES    //!< Number of properties.
    };
    
  private:
    /*! \brief Enum for the different emergent properties of molecules.
     *  \details #Emergent properties are those which need to be calculated
     *  whenever a #Property is modified. */
    enum class Emergent : size_t {
      MOLECULAR_FORMULA,    //!< Recalculate molecular formula.
      TOPOLOGICAL_BOFC,     //!< Recalculate bond orders and formal charges.
      ANGLE_PERCEPTION,     //!< Redetermine angles in molecule.
      DIHEDRAL_PERCEPTION,  //!< Redetermine dihedrals in molecule.
      NUM_EMERGENTS         //!< Number of emergent properties.
    };
    
  private:
    /*! \brief Container for storing IXAtom instances.
     *  \details A molecule takes ownership of all atoms it contains. */
    using MolAtoms = std::vector<Atom>;
    /*! \brief Container for storing IXBond instances.
     *  \details A molecule takes ownership of all bonds it contains. */
    using MolBonds = std::vector<Bond>;
    /*! \brief Container for storing IXAngle instances.
     *  \details A molecule takes ownership of all angles it contains. */
    using MolAngles = std::vector<Angle>;
    /*! \brief Container for storing IXDihedral instances.
     *  \details A molecule takes ownership of all dihedrals it contains. */
    using MolDihedrals = std::vector<Dihedral>;
    //! \brief Bitset for emergent property states.
    using EmergeSet = std::bitset<static_cast<size_>(Emergent::NUM_EMERGENTS)>;
    
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
    /*! \brief Initalisation method to always be called after construction.
     *  \details Primary reason is to work around IXMolecularGraph requiring
     *  a Molecule for construction, which can not be done until the IXMolecule
     *  has completed construction. */
    void Init();
    
    template <typename Archive>
    void Serialise(Archive& archive, const uint32_t version);
    
  public:
    /*! \brief Default constructor
     *  \details Is private to enforce that IXMolecules should only be used
     *  via the Molecule shared_ptr. */
    IXMolecule();
    
    //! \brief Destructor
    ~IXMolecule() { }
    
    /*! \brief Get the angle at position \p pos.
     *  \details Returns the angle at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of angles may
     *  change during normal operations. Angles are not re-percieved prior to
     *  retrieving the angle.
     *  \param pos the position of the angle to get.
     *  \return the angle at pos or an empty shared_ptr. */
    inline Angle GetAngle(size_ pos) const {
      return (pos < _angs.size()) ? _angs[pos] : Angle();
    }
    
    /*! \brief Get the angle at position \p pos.
     *  \details Returns the dihedral at \p pos after a range check. If \p pos
     *  is not a valid index, returns an empty shared_ptr. Positioning of
     *  dihedrals may change during normal operations. Dihedrals are not
     *  re-perceived prior to retrieving the dihedral.
     *  \param pos the position of the dihedral to get.
     *  \return the dihedral at pos or an empty shared_ptr. */
    inline Dihedral GetDihedral(size_ pos) const {
      return (pos < _dhds.size()) ? _dhds[pos] : Dihedral();
    }
    
    /*! \brief Get the angle between three atoms, if it exists.
     *  \details Searches for an angle between the three atoms and returns it.
     *  If no angle is found, the returned shared_ptr is empty. Angles are
     *  re-perceived prior to searching.
     *  \param a,b,c the atoms to get the angle between. \p b is the central
     *  atom.
     *  \return the angle between \p a and \p c centred on \p b. */
    Angle GetAngle(const Atom& a, const Atom& b, const Atom& c);
    
    /*! \brief Get the dihedral between four atoms, it if exists.
     *  \details Searches for a dihedral between the four atoms and returns it.
     *  If no dihedral is found, the returned shared_ptr is empty. Dihedrals are
     *  re-perceived prior to searching.
     *  \param a,b,c,d the atoms to get the dihedral between.
     *  \return the dihedral between \p a, \p b, \p c and \p d. */
    Dihedral GetDihedral(const Atom& a, const Atom& b,
                         const Atom& c, const Atom& d);
    
    /*! \brief Get the first angle with the given tag.
     *  \details Returns the first angle with a tag matching that given. If no
     *  such angle is found, returns an empty shared_ptr. As tags are not set on
     *  creation of angles, angles are not re-percieved prior to searching.
     *  \param tag the angle tag to search for.
     *  \return the first angle with the tag or an empty shared_ptr. */
    Angle GetAngleTag(uid_ tag) const;
    
    /*! \brief Get the first dihedral with the given tag.
     *  \details Returns the first dihedral with a tag matching that given. If
     *  no such dihedral is found, returns an empty shared_ptr. As tags are not
     *  set on creation of dihedrals, dihedrals are not re-perceived prior to
     *  searching.
     *  \param tag the dihedral tag to search for.
     *  \return the first dihedral with the tag or an empty shraed_ptr. */
    Dihedral GetDihedralTag(uid_ tag) const;
    
    /*! \brief Get the angle with the given id.
     *  \details Returns the angle with the given unique id. If no such angle is
     *  found returns an empty shared_ptr. As the unique id of any created angles
     *  will be unknown, angles are not re-percieved prior to searching.
     *  \param id the unique id of the angle to retrieve.
     *  \return the angle with the given id or an empty shared_ptr. */
    Angle GetAngleID(uid_ id) const;
    
    /*! \brief Get the dihedral with the given id.
     *  \details Returns the dihedral with the given unique id. If no such
     *  dihedral is found returns an empty shared_ptr. As the unique id of any
     *  created dihedrals will be unknown, dihedrals are not re-perceived prior
     *  to searching.
     *  \param id the unique id of the dihedral to retrieve.
     *  \return the dihedral with the given id or an empty shared_ptr. */
    Dihedral GetDihedralID(uid_ id) const;
    
    /*! \brief Get the atom at position \p pos.
     *  \details Returns that atom at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of atoms may
     *  change during normal operations.
     *  \param pos the position of the atom to get.
     *  \return the atom at pos or an empty shared_ptr. */
    inline Atom GetAtom(size_ pos) const {
      return (pos < _atms.size()) ? _atms[pos] : Atom();
    }
    
    /*! \brief Get the first atom with the given tag.
     *  \details Returns the first atom with a tag matching that given. If no
     *  such atom is found, returns an empty shared_ptr.
     *  \param tag the atom tag to search for.
     *  \return the first atom with the tag or an empty shared_ptr. */
    Atom GetAtomTag(uid_ tag) const;
    
    /*! \brief Get the atom with the given id.
     *  \details Returns the atom with the given unique id. If no such atom is
     *  found returns an empty shared_ptr.
     *  \param id the unique id of the atom to retrieve.
     *  \return the atom with the given id or an empty shared_ptr. */
    Atom GetAtomID(uid_ id) const;
    
    /*! \brief Get the bond at position \p pos.
     *  \details Returns that bond at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of bonds may
     *  change during normal operations.
     *  \param pos the position of the bond to get.
     *  \return the bond at pos or an empty shared_ptr. */
    inline Bond GetBond(size_ pos) const {
      return (pos < _bnds.size()) ? _bnds[pos] : Bond();
    }
    
    /*! \brief Get the bond between two atoms, if it exists.
     *  \details Searches for a bond between \p a and \p b and returns it. If no
     *  bond is found, returns an empty shared_ptr.
     *  \param a,b the atoms to get the bond between.
     *  \return the bond between \p a and \p b or an empty shared_ptr. */
    Bond GetBond(const Atom& a, const Atom& b) const;
    
    /*! \brief Get the first bond with the given tag.
     *  \details Returns the first bond with a tag matching that given. If no
     *  such bond is found, returns an empty shared_ptr.
     *  \param tag the bond tag to search for.
     *  \return the first bond with the tag or an empty shared_ptr. */
    Bond GetBondTag(uid_ tag) const;
    
    /*! \brief Get the bond with the given id.
     *  \details Returns the bond with the given unique id. If no such bond is
     *  found returns an empty shared_ptr.
     *  \param id the unique id of the bond to retrieve.
     *  \return the bond with the given id or an empty shared_ptr. */
    Bond GetBondID(uid_ id) const;
    
    /*! \brief Get the molecular formula of the molecule.
     *  \details Determines the molecular formula of the molecule. The elements
     *  of the formula are arranged in the order: carbon, hydrogen and then
     *  alphabetically by atomic symbol. The resultant formula is cached and
     *  recalculated only when the IXMolecule::ATOM_ELEMENTS property has been
     *  modified.
     *  \return the molecular formula of the molecule. */
    string_ GetFormula();
    
    /*! \brief Get the molecular graph for this molecule.
     *  \return the molecular graph of this molecule. */
    inline const graph::MolecularGraph& GetGraph() const { return _g; }
    
    /*! \brief Get the name of the molecule.
     *  \return the name of the molecule. */
    inline string_ GetName() const { return _name; }
    
    /*! \brief Get the molecular charge of the molecule.
     *  \return the molecular charge of the molecule. */
    inline int_ GetMolecularCharge() const { return _q; }
    
    /*! \brief Get the number of atoms in the molecule.
     *  \return the number of atoms in the molecule. */
    inline size_ NumAtoms() const { return _atms.size(); }
    
    /*! \brief Get the number of bonds in the molecule.
     *  \return the number of bonds in the molecule. */
    inline size_ NumBonds() const { return _bnds.size(); }
    
    /*! \brief Get the number of angles in the molecule.
     *  \details If the Emergent::ANGLE_PERCEPTION property is set, angles are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of angles in the molecule. */
    inline size_ NumAngles() { PerceiveAngles(); return _angs.size(); }
    
    /*! \brief Get the number of dihedrals in the molecule.
     *  \details If the CONNECTIVITY property has been modified, dihedrals are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of dihedrals in the molecule. */
    inline size_ NumDihedrals() { PerceiveDihedrals(); return _dhds.size(); }
    
    /*! \brief Set the name of the molecule.
     *  \param name the new to set. */
    inline void SetName(string_ name) { _name = name; }
    
    /*! \brief Set the molecular charge of the molecule.
     *  \details Sets the IXMolecule::ELECTRON_COUNT property as modified.
     *  \param q the new charge to set. */
    inline void SetMolecularCharge(int_ q) {
      if (q != _q) SetPropertyModified(Property::ELECTRON_COUNT);
      _q = q;
    }
    
    /*! \brief Check if the atom is owned by this molecule.
     *  \param atom the atom to check for.
     *  \return if the atom is owned by this molecule. */
    inline bool HasAtom(const Atom& atom) const {
      return atom ? atom->GetMolecule() == shared_from_this() : false;
    }
    
    /*! \brief Check if the bond is owned by this molecule.
     *  \param bond the bond to check for.
     *  \return if the bond is owned by this molecule. */
    inline bool HasBond(const Bond& bond) const {
      return bond ? bond->GetMolecule() == shared_from_this() : false;
    }
    
    /*! \brief Check if a bond between two atoms exists in this molecule.
     *  \param a,b the atoms the bond should be between.
     *  \return if there is a bond between the two atoms. */
    bool HasBond(const Atom& a, const Atom& b) const;
    
    /*! \brief Check if the angle is owned by this molecule.
     *  \param angle the angle to check for.
     *  \return if the angle is owned by this molecule. */
    inline bool HasAngle(const Angle& angle) const {
      return angle ? angle->GetMolecule() == shared_from_this() : false;
    }
    
    /*! \brief Check if an angle between three atoms exists in this molecule.
     *  \details If the Emergent::ANGLE_PERCEPTION property is set, angles are
     *  re-percieved prior to check for the existance of the angle. \p b is the
     *  central atom of the angle.
     *  \param a,b,c the atoms the angle should be between.
     *  \return if there is an angle between the three atoms. */
    bool HasAngle(const Atom& a, const Atom& b, const Atom& c);
    
    /*! \brief Check if the dihedral is owned by this molecule.
     *  \param dihedral the dihedral to check for.
     *  \return if the dihedral is owned by this molecule. */
    inline bool HasDihedral(const Dihedral& dihedral) const {
      return dihedral ? dihedral->GetMolecule() == shared_from_this() : false;
    }
    
    /*! \brief Check if a dihedral between four atoms exists in this molecule.
     *  \details If the Emergent::DIHEDRAL_PERCEPTION property is set, dihedrals
     *  are re-perceived prior to checking for the existance of the dihedral.
     *  \param a,b,c,d the atoms the dihedral should be between.
     *  \return if there is a dihedral between the four atoms. */
    bool HasDihedral(const Atom& a, const Atom& b, const Atom& c, const Atom& d);
    
    /*! \brief Create a new atom owned by the molecule.
     *  \return the new atom. */
    Atom NewAtom();
    
    /*! \brief Create a new atom of the given element for the molecule.
     *  \param element the element of the new atom.
     *  \return the new atom. */
    Atom NewAtom(Element element);
//    Atom NewAtom(Vec3 position);
    
    /*! \brief Create a new named atom owned by the molecule.
     *  \param name the name of the new atom.
     *  \return the new atom. */
    Atom NewAtom(string_ name);
    
    /*! \brief Create a new named atom of the given element for the molecule.
     *  \param name the name of the new atom.
     *  \param element the element of the new atom.
     *  \return the new atom. */
    Atom NewAtom(string_ name, Element element);
//    Atom NewAtom(string_ name, Element element, Vec3 position);
    
    /*! \brief Create a bond between two atoms.
     *  \details To create a bond, both atoms need to be owned by the molecule
     *  and there cannot be an existing bond between the atoms. If either
     *  condition is not met, the returned shared_ptr is empty.
     *  \param a,b the atoms to create a bond between.
     *  \return the new bond. */
    Bond NewBond(Atom a, Atom b);
    
    /*! \brief Remove an atom from the molecule.
     *  \details In addition to removing the atom, any bonds to the atom are
     *  also removed. A removal will not occur if the atom is not owned by the
     *  molecule it is being removed from.
     *  \param atom the atom to remove.
     *  \return if a removal occured. */
    bool RemoveAtom(Atom atom);
    
    /*! \brief Remove a bond from the molecule.
     *  \details A removal will not occur if the bond is not owned by the
     *  molecule it is being removed from. Reference to the bond is removed from
     *  both of the atoms of the bond.
     *  \param bond the bond to remove.
     *  \return if a removal occured. */
    bool RemoveBond(Bond bond);
    
    /*! \brief Remove a bond between two atoms.
     *  \param a,b the atoms to remove a bond between.
     *  \return if removal occured. */
    inline bool RemoveBond(Atom a, Atom b) { return RemoveBond(GetBond(a,b)); }
    
  private:
    /*! \brief Create an angle between three atoms.
     *  \param a,b,c the atoms to create the angle between.
     *  \return the new angle. */
    Angle NewAngle(const Atom& a, const Atom& b, const Atom& c);
    
    /*! \brief Create a dihedral between four atoms.
     *  \param a,b,c,d the atoms to create the dihedral between.
     *  \return the new dihedral. */
    Dihedral NewDihedral(const Atom& a, const Atom& b,
                         const Atom& c, const Atom& d);
    
    //! \cond
    Bond _FindBond(const Atom& a, const Atom& b) const;
    Angle _FindAngle(const Atom& a, const Atom& b, const Atom& c) const;
    Dihedral _FindDihedral(const Atom& a, const Atom& b,
                           const Atom& c, const Atom& d) const;
    //! \endcond
    
  public:
    /*! \brief Determine angles in the molecule.
     *  \details An angle is defined for each group of three connected atoms
     *  within a molecule. Subsequent calls to this method will only generate
     *  angles which have not been previously generated. Additionally, angles
     *  can only be removed by removing an atom or a bond.
     *  \return the number of angles added. */
    size_ PerceiveAngles();
    
    /*! \brief Determine dihedrals in the molecule.
     *  \details A dihedral is defined for each group of four connected atoms
     *  within a molecule. Subsequent calls to this method will only generate
     *  dihedrals which have not been previously generated. Additionally,
     *  dihedrals can only be removed by removing an atom or bond.
     *  \return the number of dihedrals added. */
    size_ PerceiveDihedrals();
//    size_t AssignElectrons();
//    bool ApplyElectronAssignment(size_t);
//    FCSCORE GetMinimumElectronAssignmentScore();
    
    /*! \brief Reserve storage space for atoms.
     *  \details Reserves storage space for a minimum of \p num IXAtom
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more atoms are added.
     *  \param num the number of IXAtoms to reserve space for. */
    inline void ReserveAtoms(size_ num) {
      if (_atms.size() < num) _atms.reserve(num);
    }
    
    /*! \brief Reserve storage space for bonds.
     *  \details Reserves storage space for a minimum of \p num IXBond
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more bonds are added.
     *  \param num the number of IXBonds to reserve space for. */
    inline void ReserveBonds(size_ num) {
      if (_bnds.size() < num ) _bnds.reserve(num);
    }
    
    /*! \brief Set that the given property has been modified.
     *  \details If a property has been modified, then any cacheable methods
     *  which require the property to calculate will be recalculated the next
     *  time they are accessed.
     *  \param prop the property that has been modified. */
    void SetPropertyModified(Property prop);
    
    /*! \brief Get iterator access to the owned atoms.
     *  \return a pair of iterators indicating the begining and end of the
     *  owned atoms. */
    inline std::pair<MolAtomIter, MolAtomIter> GetAtoms() const {
      return std::make_pair(_atms.begin(), _atms.end());
    }
    
    /*! \brief Get iterator access to the owned bonds.
     *  \return a pair of iterators indication the beginning and end of the
     *  owned bonds. */
    inline std::pair<MolBondIter, MolBondIter> GetBonds() const {
      return std::make_pair(_bnds.begin(), _bnds.end());
    }
    
    /*! \brief Get iterator access to the owned angles.
     *  \details Angles are re-percieved prior to returning the iterators.
     *  \return a pair of iterators for the beginning and end of the angles. */
    inline std::pair<MolAngleIter, MolAngleIter> GetAngles() {
      PerceiveAngles();
      return std::make_pair(_angs.begin(), _angs.end());
    }
    
    /*! \brief Get iterator access to the owned dihedrals.
     *  \details Dihedrals are re-perceived prior to returning the iterators.
     *  \return a pair of iterators indication the beginning and end of the
     *  owned dihedrals. */
    inline std::pair<MolDihedralIter, MolDihedralIter> GetDihedrals() {
      PerceiveDihedrals();
      return std::make_pair(_dhds.begin(), _dhds.end());
    }
    
  private:
    //! Name of the molecule
    string_ _name;
    //! Molecular charge of the molecule
    int_ _q;
    //! Atoms owned by molecule
    MolAtoms _atms;
    //! Bonds owned by molecule
    MolBonds _bnds;
    //! Angles owned by molecule
    MolAngles _angs;
    //! Dihedrals owned by molecule
    MolDihedrals _dhds;
    //! Mask for modified emergent properties
    EmergeSet _emerge;
    //! Molecular graph of molecule
    graph::MolecularGraph _g;
    //! Cached molecular formula string
    string_ _formula_cache;
  };
  
  /*! \brief Creates a new molecule.
   *  \details Strictly enforces usage of the IXMolecule class through use of
   *  the Molecule shared_ptr.
   *  \return a new, blank Molecule. */
  inline Molecule CreateMolecule() {
    Molecule m(std::make_shared<IXMolecule>()); m->Init(); return m; }
  
  //! \brief Type for the modifiable properties of an IXMolecule
  using MolProperty = indigox::IXMolecule::Property;
  
}

#endif /* INDIGOX_MOLECULE_HPP */
