/*! \file molecule.hpp */
#include <bitset>
#include <map>
#include <unordered_map>
#include <vector>

#include "../utils/counter.hpp"
#include "../utils/numerics.hpp"

#ifndef INDIGOX_MOLECULE_HPP
#define INDIGOX_MOLECULE_HPP

namespace indigox {
  
  class IXAtom;
  class IXBond;
  class IXAngle;
  class IXDihedral;
  class IXMolecule;
  namespace test { class IXMolecule; }
  class IXElement;
  namespace graph {
    class IXMolecularGraph;
    typedef std::shared_ptr<IXMolecularGraph> MolecularGraph;
  }
  
  typedef std::shared_ptr<IXAtom> Atom;
  typedef std::shared_ptr<IXBond> Bond;
  typedef std::shared_ptr<IXAngle> Angle;
  typedef std::shared_ptr<IXDihedral> Dihedral;
  //! \brief shared_ptr for normal use of the IXMolecule class.
  typedef std::shared_ptr<IXMolecule> Molecule;
  typedef std::shared_ptr<IXElement> Element;
  
  typedef std::weak_ptr<IXAtom> _Atom;
  typedef std::weak_ptr<IXBond> _Bond;
  typedef std::weak_ptr<IXAngle> _Angle;
  typedef std::weak_ptr<IXDihedral> _Dihedral;
  /*! \brief weak_ptr for non-ownership reference to the IXMolecule class.
   *  \details Intended for internal use only. */
  typedef std::weak_ptr<IXMolecule> _Molecule;
  typedef std::weak_ptr<IXElement> _Element;
  
  class IXMolecule : public utils::IXCountableObject<IXMolecule>,
  public std::enable_shared_from_this<IXMolecule> {
    //! \brief Friendship allows generation of molecules.
    friend Molecule CreateMolecule();
    //! \brief Friendship allows IXMolecule to be tested.
    friend class indigox::test::IXMolecule;
    
  private:
    /*! \brief Container for storing IXAtom instances.
     *  \details A molecule takes ownership of all atoms it contains. */
    typedef std::vector<Atom> MolAtoms;
    /*! \brief Container for storing IXBond instances.
     *  \details A molecule takes ownership of all bonds it contains. */
    typedef std::vector<Bond> MolBonds;
    /*! \brief Container for storing IXAngle instances.
     *  \details A molecule takes ownership of all angles it contains. */
    typedef std::vector<Angle> MolAngles;
    /*! \brief Container for storing IXDihedral instances.
     *  \details A molecule takes ownership of all dihedrals it contains. */
    typedef std::vector<Dihedral> MolDihedrals;
    
  public:   // Public iterator typedefs for easier external usage
    //! \brief Iterator over owned IXAtom instances.
    typedef MolAtoms::iterator MolAtomIter;
    //! \brief Iterator over owned IXBond instances.
    typedef MolBonds::iterator MolBondIter;
    //! \brief Iterator over owned IXAngle instances.
    typedef MolAngles::iterator MolAngleIter;
    //! \brief Iterator over owned IXDihedral instances.
    typedef MolDihedrals::iterator MolDihedralIter;
    
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
    /*! \brief Default constructor
     *  \details Is private to enforce that IXMolecules should only be used
     *  via the Molecule shared_ptr. */
    IXMolecule();
    
    /*! \brief Initalisation method to always be called after construction.
     *  \details Primary reason is to work around IXMolecularGraph requiring
     *  a Molecule for construction, which can not be done until the IXMolecule
     *  has completed construction. */
    void Init();
    
  public:
    //! \brief Destructor
    ~IXMolecule() { }
    
    /*! \brief Get the atom at position \p pos.
     *  \details Returns that atom at \p pos after a range check. If \p pos is
     *  not a valid index, returns an empty shared_ptr. Positioning of atoms may
     *  change during normal operations.
     *  \param pos the position of the atom to get.
     *  \return the atom at pos or an empty shared_ptr. */
    inline Atom GetAtom(size_ pos) const {
      return (pos < _atoms.size()) ? _atoms[pos] : Atom();
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
    Bond GetBond(size_ pos) const {
      return (pos < _bonds.size()) ? _bonds[pos] : Bond();
    }
    
    /*! \brief Get the bond between two atoms, if it exists.
     *  \details Searches for a bond between \p a and \p b and returns it. If no
     *  bond is found, returns an empty shared_ptr.
     *  \param a,b the atoms to get the bond between.
     *  \return the bond between \p a and \p b or an empty shared_ptr. */
    Bond GetBond(Atom a, Atom b) const;
    
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
    graph::MolecularGraph GetGraph() const { return _g; }
    
    /*! \brief Get the name of the molecule.
     *  \return the name of the molecule. */
    string_ GetName() const { return _name; }
    
    /*! \brief Get the molecular charge of the molecule.
     *  \return the molecular charge of the molecule. */
    int GetMolecularCharge() const { return _q; }
    
    /*! \brief Get the number of atoms in the molecule.
     *  \return the number of atoms in the molecule. */
    size_ NumAtoms() const { return _atoms.size(); }
    
    /*! \brief Get the number of bonds in the molecule.
     *  \return the number of bonds in the molecule. */
    size_ NumBonds() const { return _bonds.size(); }
    
    /*! \brief Get the number of angles in the molecule.
     *  \details If the CONNECTIVITY property has been modified, angles are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of angles in the molecule. */
    size_ NumAngles();
    
    /*! \brief Get the number of dihedrals in the molecule.
     *  \details If the CONNECTIVITY property has been modified, dihedrals are
     *  re-calculated prior to determining how many there are in the molecule.
     *  \return the number of dihedrals in the molecule. */
    size_ NumDihedrals();
    
    /*! \brief Set the name of the molecule.
     *  \param name the new to set. */
    void SetName(string_ name) { _name = name; }
    
    /*! \brief Set the molecular charge of the molecule.
     *  \details Sets the IXMolecule::ELECTRON_COUNT property as modified.
     *  \param q the new charge to set. */
    void SetMolecularCharge(int q);
    
    /*! \brief Check if the atom is owned by this molecule.
     *  \param atom the atom to check for.
     *  \return if the atom is owned by this molecule. */
    bool HasAtom(Atom atom) const;
    
    /*! \brief Check if the bond is owned by this molecule.
     *  \param bond the bond to check for.
     *  \return if the bond is owned by this molecule. */
    bool HasBond(Bond bond) const;
    
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
    
//    size_t AssignElectrons();
//    bool ApplyElectronAssignment(size_t);
//    FCSCORE GetMinimumElectronAssignmentScore();
    
    /*! \brief Reserve storage space for atoms.
     *  \details Reserves storage space for a minimum of \p num IXAtom
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more atoms are added.
     *  \param num the number of IXAtoms to reserve space for. */
    void ReserveAtoms(size_ num) {
      if (_atoms.size() < num) _atoms.reserve(num);
    }
    
    /*! \brief Reserve storage space for bonds.
     *  \details Reserves storage space for a minimum of \p num IXBond
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more bonds are added.
     *  \param num the number of IXBonds to reserve space for. */
    void ReserveBonds(size_ num) {
      if (_bonds.size() < num ) _bonds.reserve(num);
    }
    
    /*! \brief Reserve storage space for angles.
     *  \details Reserves storage space for a minimum of \p num IXAngle
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more angles are added.
     *  \param num the number of IXAngles to reserve space for. */
    void ReserveAngles(size_ num) {
      if (_angles.size() < num) _angles.reserve(num);
    }
    
    /*! \brief Reserve storage space for dihedrals.
     *  \details Reserves storage space for a minimum of \p num IXDihedral
     *  instances. This is more efficient when building large molecules as the
     *  vector will not need to grow as more dihedrals are added.
     *  \param num the number of IXDihedrals to reserve space for. */
    void ReserveDihedrals(size_ num) {
      if (_dihedrals.size() < num) _dihedrals.reserve(num);
    }
    
    /*! \brief Set that the given property has been modified.
     *  \details If a property has been modified, then any cacheable methods
     *  which require the property to calculate will be recalculated the next
     *  time they are accessed.
     *  \param prop the property that has been modified. */
    void SetPropertyModified(Property prop);
    
  private:
    //! Name of the molecule
    string_ _name;
    //! Molecular charge of the molecule
    int_ _q;
    //! Atoms owned by molecule
    MolAtoms _atoms;
    //! Bonds owned by molecule
    MolBonds _bonds;
    //! Angles owned by molecule
    MolAngles _angles;
    //! Dihedrals owned by molecule
    MolDihedrals _dihedrals;
    //! Mask for modified emergent properties
    std::bitset<static_cast<size_>(Emergent::NUM_EMERGENTS)> _emerge;
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
    Molecule m(new indigox::IXMolecule()); m->Init(); return m; }
  
}

#endif /* INDIGOX_MOLECULE_HPP */
