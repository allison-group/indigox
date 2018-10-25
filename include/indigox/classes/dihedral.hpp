/*! \file dihedral.hpp */
#include <array>
#include <cstdint>
#include <memory>
#include <string>

#include "../utils/counter.hpp"
#include "../utils/fwd_declares.hpp"
#include "../utils/quad.hpp"

#ifndef INDIGOX_CLASSES_DIHEDRAL_HPP
#define INDIGOX_CLASSES_DIHEDRAL_HPP

namespace indigox {
  class Dihedral
  : public utils::IXCountableObject<Dihedral>,
  public std::enable_shared_from_this<Dihedral> {
    //! \brief Friendship allows Molecule to create new dihedrals.
    friend class indigox::Molecule;
    //! \brief Friendship allows Dihedral to be tested.
    friend struct indigox::test::TestDihedral;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
  private:
    //! \brief Container for storing IXAtom references assigned to an IXDihedral.
    using DihedAtoms = std::array<wAtom, 4>;
    
  public:
    //! \brief Type for storing dihedral parameters
    using DihedTypes = std::vector<wFFDihedral>;
    //! \brief Iterator over IXAtom references stored on an IXDihedral.
    using DihedAtomIter = DihedAtoms::const_iterator;
    
  private:
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<Dihedral>& construct,
                                   const uint32_t version);
    void Clear();
  public:
    Dihedral() = delete;  // no default constructor
    
    /*! \brief Normal constructor.
     *  \details Creates a dihedral between four atoms, linking it to the given
     *  Molecule.
     *  \param a,b,c,d the four atoms to construct a dihedral between.
     *  \param m the molecule to assign the angle to. */
    Dihedral(Atom& a, Atom& b, Atom& c, Atom& d, Molecule& m);
    
    //! \brief Destructor
    ~Dihedral() { }
    
    /*! \brief Tag of the dihedral.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the dihedral, use IXDihedral::GetUniqueID().
     *  \return the tag assigned to the dihedral. */
    inline uint32_t GetTag() const { return _tag; }
    
    /*! \brief Molecule this dihedral is associated with.
     *  \details The returned shared_ptr is empty if the dihedral is not
     *  assigned to a valid molecule.
     *  \return the molecule associated with this dihedral. */
    inline Molecule& GetMolecule() const { return *_mol.lock(); }
    
    /*! \brief Get the atoms of the dihedral.
     *  \return quad of the atoms of the dihedral. */
    inline stdx::quad<Atom&, Atom&, Atom&, Atom&> GetAtoms() const {
      return {*_atms[0].lock(), *_atms[1].lock(),
              *_atms[2].lock(), *_atms[3].lock()};
    }
    
    /*! \brief Number of atoms this dihedral contains.
     *  \return 4. */
    size_t NumAtoms() const { return _atms.size(); }
    
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
    std::string ToString() const;
    
    /*! \brief Set the tag of this dihedral.
     *  \details The tag of a dihedral should not be considered stable. Use with
     *  caution.
     *  \param tag the tag to set. */
    inline void SetTag(uint32_t tag) { _tag = tag; }
    
    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the dihedral in the container of
     *  dihedrals of the molecule it is a part of. If the molecule is dead, the
     *  index returned is the tag of the dihedral.
     *  \return the index of the dihedral. */
    size_t GetIndex() const;
    
    /*! \brief Get the type of the dihedral.
     *  \return the type of the dihedral. */
    FFDihedral& GetType(size_t pos) const;
    size_t NumTypes() const { return _types.size(); }
    const DihedTypes& GetTypes() const { return _types; }
    
    
    /*! \brief Set the type of the dihedral.
     *  \param type the type of dihedral to set. */
    void AddType(FFDihedral& type);
    void RemoveType(FFDihedral& type);
    
  private:
    //! The molecule this dihedral is assigned to.
    wMolecule _mol;
    //! Tag (unstable)
    uint32_t _tag;
    //! \brief Atoms which make up the dihedral.
    DihedAtoms _atms;
    //! \brief Type of dihedral
    DihedTypes _types;
  };
  
  std::ostream& operator<<(std::ostream& os, const Dihedral& dhd);
}

#endif /* INDIGOX_CLASSES_DIHEDRAL_HPP */
