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
  class Dihedral {
    //! \brief Friendship allows Molecule to create new dihedrals.
    friend class indigox::Molecule;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
  public:
    //! \brief Container for storing IXAtom references assigned to an IXDihedral.
    using DihedralAtoms = std::array<Atom, 4>;
    //! \brief Type for storing dihedral parameters
    using DihedralTypes = std::vector<FFDihedral>;
    
  private:
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(Dihedral);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(Dihedral, dhd);
    
  private:
    /*! \brief Normal constructor.
     *  \details Creates a dihedral between four atoms, linking it to the given
     *  Molecule.
     *  \param a,b,c,d the four atoms to construct a dihedral between.
     *  \param m the molecule to assign the angle to. */
    Dihedral(const Atom& a, const Atom& b, const Atom& c, const Atom& d,
             const Molecule& m);
    
  public:
    /*! \brief Tag of the dihedral.
     *  \details This value may be modified without warning. Use with caution.
     *  For a constant identifier to the dihedral, use IXDihedral::GetUniqueID().
     *  \return the tag assigned to the dihedral. */
    int32_t GetTag() const;
    
    /*! \brief Molecule this dihedral is associated with.
     *  \details The returned shared_ptr is empty if the dihedral is not
     *  assigned to a valid molecule.
     *  \return the molecule associated with this dihedral. */
    const Molecule& GetMolecule() const;
    
    /*! \brief Get the atoms of the dihedral.
     *  \return quad of the atoms of the dihedral. */
    const DihedralAtoms& GetAtoms() const;
    
    /*! \brief Number of atoms this dihedral contains.
     *  \return 4. */
    constexpr int64_t NumAtoms() const { return 4; }
    
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
    void SetTag(int32_t tag);
    
    /*! \brief Get the index from the molecule.
     *  \details Calculates the index of the dihedral in the container of
     *  dihedrals of the molecule it is a part of. If the molecule is dead, the
     *  index returned is the tag of the dihedral.
     *  \return the index of the dihedral. */
    int64_t GetIndex() const;
    
    /*! \brief Get the type of the dihedral.
     *  \return the type of the dihedral. */
    int64_t NumTypes() const;
    const DihedralTypes& GetTypes() const;
    bool HasType() const;
    void SetTypes(const DihedralTypes& types)
//    {
//      _types.clear();
//      _types.assign(types.begin(), types.end());
//    }
    /*! \brief Set the type of the dihedral.
     *  \param type the type of dihedral to set. */
    void AddType(const FFDihedral& type);
    void RemoveType(const FFDihedral& type);
    int32_t GetPriority() const;
    
  private:
    struct Impl;
    std::shared_ptr<Impl> m_data;
    
    //! The molecule this dihedral is assigned to.
    wMolecule _mol;
    //! Tag (unstable)
    uint32_t _tag;
    //! \brief Atoms which make up the dihedral.
    DihedAtoms _atms;
    //! \brief Type of dihedral
    DihedTypes _types;
    //! \brief Priority of dihedral. Only set by PerceiveDihedrals
    int32_t _priority;
  };
  
  std::ostream& operator<<(std::ostream& os, const Dihedral& dhd);
}

#endif /* INDIGOX_CLASSES_DIHEDRAL_HPP */
