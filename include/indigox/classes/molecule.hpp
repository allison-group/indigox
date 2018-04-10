/** @file molecule.hpp
 *  @brief Declaration of Molecule class
 *  @author Ivan Welsh
 *  @date 6 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */
#include <map>
#include <unordered_map>
#include <vector>

#include "../algorithm/electron_optimisation.hpp"

#include "molecular_graph.hpp"
#include "../utils/counter.hpp"

#ifndef INDIGOX_MOLECULE_HPP
#define INDIGOX_MOLECULE_HPP

namespace indigox {
  
  // Related typedefs
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
  typedef std::vector<Atom> MolAtoms;
  typedef std::vector<Bond> MolBonds;
  typedef MolAtoms::iterator MolAtomIterator;
  typedef MolBonds::iterator MolBondIterator;
  typedef MolAtoms::const_iterator const_MolAtomIterator;
  typedef MolBonds::const_iterator const_MolBondIterator;
  
  /** @class Molecule molecule.hpp classes/molecule.hpp
   *  @brief Class for storing and manipulating molecules.
   */
  class IXMolecule : public utils::CountableObject<IXMolecule>,
  public std::enable_shared_from_this<IXMolecule> {
    
  public:
    /// @name Initialisation methods
    
    /// @brief Default constructor
    IXMolecule();
    
    /// @brief Constructor for pre-named molecule
    IXMolecule(std::string name);
    
    /// @brief Destructor
    ~IXMolecule();
    
    /// @name Data retrival methods
    
    /// @returns a shared pointer to an atom with the given index.
    Atom GetAtomIndex(uid_t);
    Atom GetAtomUniqueID(uid_t);
    
    /// @returns a shared pointer to a bond with the given index.
    Bond GetBondIndex(uid_t);
    Bond GetBondUniqueID(uid_t);
    
    /// @returns a shared pointer to the bond between two atoms, if it exists.
    Bond GetBond(Atom, Atom) const;
    
    /// @returns the chemical formula of the molecule.
    std::string GetFormula() const;
    
    /// @returns a shared pointer to a graph representation of the molecule
    MolecularGraph GetMolecularGraph();
    
    /// @returns the name of the molecule.
    std::string GetName() const;
    
    /// @returns the molecular charge of the molecule
    int GetTotalCharge() const;
    
    bool IsModified() const;
    
    /// @returns the number of atoms in the molecule.
    size_t NumAtoms() const;
    
    /// @returns the number of bonds in the molecule.
    size_t NumBonds() const;
    
    /// @name Data modification methods
    
    /// @brief Sets the name of the molecule
    void SetName(std::string);
    
    /// @brief Sets the molecular charge of the molecule.
    void SetTotalCharge(int);
    
    /// @brief Resets the indices of contained atoms and bonds.
    void ResetIndices();
    
    /// @name Molecule modification methods
    
    /// @brief Adds a new atom to the molecule.
    /// @returns a shared pointer to the new atom.
    Atom NewAtom();
    Atom NewAtom(Element);
    Atom NewAtom(uid_t, Element);
    
    /// @brief Adds a new bond between two atoms to the molecule.
    /// @returns a shared pointer to the new bond.
    Bond NewBond(Atom, Atom);
    
    /// @brief Removes a given atom from the molecule.
    void RemoveAtom(Atom);
    
    /// @brief Removes a given bond from the molecule.
    void RemoveBond(Bond);
    
    size_t AssignElectrons();
    bool ApplyElectronAssignment(size_t);
    FCSCORE GetMinimumElectronAssignmentScore();
    
    /// @name Iterator methods
    Atom Begin(MolAtomIterator&);
    Atom Next(MolAtomIterator&);
    
    MolAtomIterator BeginAtom();
    MolAtomIterator EndAtom();
    MolBondIterator BeginBond();
    MolBondIterator EndBond();
    
    
  private:
    std::string name_;
    int q_ = 0;
    MolAtoms atoms_;  // own atoms that are part of me
    MolBonds bonds_;  // own bonds that are part of me
    std::unordered_map<uid_t, std::weak_ptr<IXAtom>> idx_to_atom_;
    std::unordered_map<uid_t, std::weak_ptr<IXBond>> idx_to_bond_;
    std::map<Atom, MolVertex> atom_to_vertex_;
    std::map<Bond, MolEdge> bond_to_edge_;
    MolecularGraph graph_;
    std::unique_ptr<ElectronOpt> elnopt_;
    bool modified_ = true;
  };
  
}

#endif /* INDIGOX_MOLECULE_HPP */
