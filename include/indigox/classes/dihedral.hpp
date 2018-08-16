/*! \file dihedral.hpp */
#include <array>
#include <memory>

#include "../utils/counter.hpp"
#include "../utils/numerics.hpp"
#include "../utils/quad.hpp"

#ifndef INDIGOX_CLASSES_DIHEDRAL_HPP
#define INDIGOX_CLASSES_DIHEDRAL_HPP

namespace indigox {
  class IXAtom;
  class IXDihedral;
  class IXMolecule;
  class IXFFDihedral;
  namespace test { struct TestDihedral; }
  
  using Atom = std::shared_ptr<IXAtom>;
  //! \brief shared_ptr for normal use of the IXDihedral class.
  using Dihedral = std::shared_ptr<IXDihedral>;
  using Molecule = std::shared_ptr<IXMolecule>;
  using FFDihedral = std::shared_ptr<IXFFDihedral>;
  
  using _Atom = std::weak_ptr<IXAtom>;
  /*! \brief weak_ptr for non-ownership reference to the IXDihedral class.
   *  \details Intended for internal use only. */
  using _Dihedral = std::weak_ptr<IXDihedral>;
  using _Molecule = std::weak_ptr<IXMolecule>;
  
  class IXDihedral
  : public utils::IXCountableObject<IXDihedral>,
  public std::enable_shared_from_this<IXDihedral> {
    //! \brief Friendship allows IXMolecule to create new dihedrals.
    friend class indigox::IXMolecule;
    //! \brief Friendship allows IXDihedral to be tested.
    friend struct indigox::test::TestDihedral;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
  private:
    //! \brief Container for storing IXAtom references assigned to an IXDihedral.
    using DihedAtoms = std::array<_Atom, 4>;
    
  public:
    //! \brief Iterator over IXAtom references stored on an IXDihedral.
    using DihedAtomIter = DihedAtoms::const_iterator;
    
  private:
    /*! \brief Normal constructor.
     *  \details Creates a dihedral between four atoms, linking it to the given
     *  Molecule.
     *  \param a,b,c,d the four atoms to construct a dihedral between.
     *  \param m the molecule to assign the angle to. */
    IXDihedral(Atom a, Atom b, Atom c, Atom d, Molecule m);
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXDihedral>& construct,
                                   const uint32_t version);
    
  public:
    IXDihedral() = delete;  // no default constructor
    
    //! \brief Destructor
    ~IXDihedral() { }
    
    /*! \brief Tag of the dihedral.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the dihedral, use IXDihedral::GetUniqueID().
     *  \return the tag assigned to the dihedral. */
    inline uid_ GetTag() const { return _tag; }
    
    /*! \brief Molecule this dihedral is associated with.
     *  \details The returned shared_ptr is empty if the dihedral is not
     *  assigned to a valid molecule.
     *  \return the molecule associated with this dihedral. */
    inline Molecule GetMolecule() const { return _mol.lock(); }
    
    /*! \brief Get the atoms of the dihedral.
     *  \return quad of the atoms of the dihedral. */
    inline stdx::quad<Atom, Atom, Atom, Atom> GetAtoms() const {
      return stdx::make_quad(_atms[0].lock(), _atms[1].lock(),
                             _atms[2].lock(), _atms[3].lock());
    }
    
    /*! \brief Number of atoms this dihedral contains.
     *  \return 4. */
    size_ NumAtoms() const { return _atms.size(); }
    
    /*! \brief Switch the order of atoms in the dihedral.
     *  \details The two outer atoms of the dihedral are swapped, as are the two
     *  inner atoms. */
    inline void SwapOrder() {
      std::swap(_atms[0], _atms[3]);
      std::swap(_atms[1], _atms[2]);
    }
    
    /*! \brief String representation of the dihedral.
     *  \details If any of the atoms of the dihedral are not valid, the returned
     *  string is Dihedral(MALFORMED), otherwise it is of the form:
     *  Dihedral(NAME_A, NAME_B, NAME_C, NAME_D).
     *  \return a string representation of the dihedral. */
    string_ ToString() const;
    
    /*! \brief Set the tag of this dihedral.
     *  \details The tag of a dihedral should not be considered stable. Use with
     *  caution.
     *  \param tag the tag to set. */
    inline void SetTag(uid_ tag) { _tag = tag; }
    
    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the dihedral in the container of
     *  dihedrals of the molecule it is a part of. If the molecule is dead, the
     *  index returned is the tag of the dihedral.
     *  \return the index of the dihedral. */
    size_ GetIndex() const;
    
    /*! \brief Get the type of the dihedral.
     *  \return the type of the dihedral. */
    FFDihedral GetType() const { return _type; }
    
    /*! \brief Set the type of the dihedral.
     *  \param type the type of dihedral to set. */
    void SetType(FFDihedral type) { _type = type; }
    
  private:
    /*! \brief Clear all informations.
     *  \details Erases all information stored on the angle, and resets
     *  everything back to a just created state. */
    void Clear();
    
  private:
    //! The molecule this dihedral is assigned to.
    _Molecule _mol;
    //! Tag (unstable)
    uid_ _tag;
    //! \brief Atoms which make up the dihedral.
    DihedAtoms _atms;
    //! \brief Type of dihedral
    FFDihedral _type;
  };
  
  std::ostream& operator<<(std::ostream& os, const IXDihedral& dhd);
  inline std::ostream& operator<<(std::ostream& os, const Dihedral& dhd) {
    return dhd ? os << *dhd : os;
  }
}

#endif /* INDIGOX_CLASSES_DIHEDRAL_HPP */
