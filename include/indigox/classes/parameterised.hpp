#include <array>
#include <vector>

#include <EASTL/vector_map.h>

#include "../utils/fwd_declares.hpp"
#include "../utils/numerics.hpp"
#include "../utils/triple.hpp"
#include "../utils/quad.hpp"

#ifndef INDIGOX_CLASSES_PARAMETERISED_HPP
#define INDIGOX_CLASSES_PARAMETERISED_HPP

namespace indigox {
  
  class ParamAtom {
  public:
    //! \brief Type giving counts of each mapped atom type
    using TypeCounts = eastl::vector_map<FFAtom, size_t>;
    //! \brief Type giving all mapped charges
    using MappedCharge = std::vector<double>;
    //! \brief Type giving all the mapped atoms
    using MappedAtoms = std::vector<wAtom>;
    friend class ParamMolecule;
    
  public:
    ParamAtom();
    ParamAtom(const ParamAtom& atm);
    ParamAtom(ParamAtom&& atm);
    ParamAtom& operator=(const ParamAtom& atm);
    ParamAtom& operator=(ParamAtom&& atm);
    
  private:
    /*! \brief Normal constructor
     *  \param atm the Atom this is parameterising. */
    ParamAtom(Atom& atm);
    
  public:
    /*! \brief Obtain details from mapped atom
     *  \param mapped the atom matched. */
    void MappedWith(Atom& mapped);
    
    /*! \brief Apply the parameterisation.
     *  \details Applies the parameteristion. Doing so sets the partial charge
     *  on the atom to the mean of the charge, and the atom type to the most
     *  common type. Additionally, doing so will mean that no more atoms can be
     *  mapped with. If the \p self_consistent flag is true, there is a
     *  requirement that all the mapped atoms have the same parameters. As such,
     *  if there arises a situtation where this requirement is not met, an
     *  exception will be thrown.
     *  \param self_consistent if the parameterisation needs to be self
     *  consistent. */
    bool ApplyParameterisation(bool self_consistent);
    
  public:
    /*! \brief Calculates the number of mapped atoms
     *  \return the number of atoms which have currently been mapped to me. */
    int64_t NumSourceAtoms() const;
    
    /*! \brief Get the atom this parameterises.
     *  \return the atom this parameterises. */
    Atom& GetAtom() const;
    
    /*! \brief Get the mean of the charges
     *  \return the mean of the charges. */
    double MeanCharge() const;
    
    /*! \brief Get the median of the charges.
     *  \return the median of the charges. */
    double MeadianCharge();
    
    /*! \brief Get the standard deviation of the charges.
     *  \return the standard deviation of the charges. */
    double StandardDeviationCharge() const;
    
    /*! \brief Get the mode of mapped atom types.
     *  \return the most commonly mapped atom type. */
    FFAtom GetMostCommonType() const;
    
    const TypeCounts& GetMappedTypeCounts() const;
    const MappedCharge& GetMappedCharges() const;
    
  public:
    bool operator==(const ParamAtom& atm) const;
    bool operator!=(const ParamAtom& atm) const;
    bool operator<(const ParamAtom& atm) const;
    bool operator>(const ParamAtom& atm) const;
    bool operator<=(const ParamAtom& atm) const;
    bool operator>=(const ParamAtom& atm) const;
    operator bool() const;
    
  private:
    struct ParamAtomImpl;
    std::shared_ptr<ParamAtomImpl> m_patmdat;
  };
  std::ostream& operator<<(std::ostream& os, const ParamAtom& atm);
  
  class ParamBond {
  public:
    //! \brief Type giving counts of each mapped bond type
    using TypeCounts = eastl::vector_map<FFBond, size_t>;
    //! \brief Type giving all the mapped bonds
    using MappedBonds = std::vector<wBond>;
    //! \brief Type giving the atoms this parameterised bond is between
    using BondAtoms = std::pair<wAtom, wAtom>;
    friend class ParamMolecule;
    
  public:
    ParamBond();
    ParamBond(const ParamBond& bnd);
    ParamBond(ParamBond&& bnd);
    ParamBond& operator=(const ParamBond& bnd);
    ParamBond& operator=(ParamBond&& bnd);
    
  private:
    /*! \brief Normal constructor
     *  \param a,b the atoms which mark the bond this is parameterising. */
    ParamBond(std::pair<Atom&, Atom&> atms, Bond& bnd);
    
  public:
    /*! \brief Obtain details from mapped bonds.
     *  \param mapped the bond matched. */
    void MappedWith(Bond& mapped);
    
    /*! \brief Apply the parameterisation.
     *  \details Applies the parameteristion. Doing so sets the bond type to the
     *  most common type. Additionally, doing so will mean that no more bonds
     *  can be mapped with. If the \p self_consistent flag is true, there is a
     *  requirement that all the mapped bonds have the same parameters. As such,
     *  if there arises a situtation where this requirement is not met, an
     *  exception will be thrown.
     *  \param self_consistent if the parameterisation needs to be self
     *  consistent. */
    bool ApplyParameterisation(bool self_consistent);
    
  public:
    /*! \brief Calculates the number of mapped bonds.
     *  \return the number of bonds which have currently been mapped. */
    int64_t NumSourceBonds() const;
    
    /*! \brief Get the atoms this parameterises the bond between.
     *  \return the two atoms defining the bond this parameterises. */
    std::pair<Atom&, Atom&> GetAtoms() const;
    
    /*! \brief Get the bond that is parameterised.
     *  \return the parameterised bond. */
    Bond& GetBond() const;
    
    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    FFBond GetMostCommonType() const;
    
    const TypeCounts& GetMappedTypeCounts() const;
    
  public:
    bool operator==(const ParamBond& bnd) const;
    bool operator!=(const ParamBond& bnd) const;
    bool operator<(const  ParamBond& bnd) const;
    bool operator>(const  ParamBond& bnd) const;
    bool operator<=(const ParamBond& bnd) const;
    bool operator>=(const ParamBond& bnd) const;
    operator bool() const;
    
  private:
    struct ParamBondImpl;
    std::shared_ptr<ParamBondImpl> m_pbnddat;
  };
  std::ostream& operator<<(std::ostream& os, const ParamBond& bnd);
  
  class ParamAngle {
  public:
    //! \brief Type giving counts of each mapped angle type
    using TypeCounts = eastl::vector_map<FFAngle, size_t>;
    //! \brief Type giving all the mapped angles
    using MappedAngles = std::vector<wAngle>;
    //! \brief Type giving the atoms this parameterised angle is between
    using AngleAtoms = stdx::triple<wAtom>;
    friend class ParamMolecule;
    
  public:
    ParamAngle();
    ParamAngle(const ParamAngle& ang);
    ParamAngle(ParamAngle&& ang);
    ParamAngle& operator=(const ParamAngle& ang);
    ParamAngle& operator=(ParamAngle& ang);
    
  private:
    /*! \brief Normal constructor.
     *  \param a,b,c the atoms which mark the parameterised angle. */
    ParamAngle(stdx::triple<Atom&> atms, Angle& ang);
    
  public:
    /*! \brief Obtain details from mapped angles.
     *  \param mapped the angle matched. */
    void MappedWith(Angle& mapped);
    
    /*! \brief Apply the parameterisation.
     *  \details Applies the parameteristion. Doing so sets the angle type to
     *  the most common type. Additionally, doing so will mean that no more
     *  angles can be mapped with. If the \p self_consistent flag is true, there
     *  is a requirement that all the mapped angles have the same parameters. As
     *  such, if there arises a situtation where this requirement is not met, an
     *  exception will be thrown.
     *  \param self_consistent if the parameterisation needs to be self
     *  consistent. */
    bool ApplyParameterisation(bool self_consistent);
    
  public:
    /*! \brief Calculate the number of mapped angles.
     *  \return the number of angles which have currently been mapped. */
    int64_t NumSourceAngles() const;
    
    /*! \brief Get the atoms this parameterises the angle between.
     *  \return the three atoms defining the angle this parameterises. */
    stdx::triple<Atom&> GetAtoms() const;
    
    /*! \brief Get the angle that is parameterised.
     *  \return the parameterised bond. */
    Angle& GetAngle() const;
    
    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    FFAngle GetMostCommonType() const;
    
    const TypeCounts& GetMappedTypeCounts() const;
    
  public:
    bool operator==(const ParamAngle& ang) const;
    bool operator!=(const ParamAngle& ang) const;
    bool operator<(const  ParamAngle& ang) const;
    bool operator>(const  ParamAngle& ang) const;
    bool operator<=(const ParamAngle& ang) const;
    bool operator>=(const ParamAngle& ang) const;
    operator bool() const;
    
  private:
    struct ParamAngleImpl;
    std::shared_ptr<ParamAngleImpl> m_pangdat;
  };
  std::ostream& operator<<(std::ostream& os, const ParamAngle& ang);
  
  class ParamDihedral {
  public:
    //! \brief Type giving counts of each mapped dihedral type
    using TypeGroup = std::vector<FFDihedral>;
    using TypeCounts = eastl::vector_map<TypeGroup, size_t>;
    //! \brief Type giving all the mapped dihedrals
    using MappedDihedrals = std::vector<wDihedral>;
    //! \brief Type giving the atoms this parameterised dihedral is between
    using DihedralAtoms = stdx::quad<wAtom>;
    
    friend class ParamMolecule;
    
  public:
    ParamDihedral();
    ParamDihedral(const ParamDihedral& dhd);
    ParamDihedral(ParamDihedral&& dhd);
    ParamDihedral& operator=(const ParamDihedral& dhd);
    ParamDihedral& operator=(ParamDihedral&& dhd);
    
  private:
    /*! \brief Normal constructor.
     *  \param a,b,c,d the atoms which mark the parameterised dihedral. */
    ParamDihedral(stdx::quad<Atom&> atms, Dihedral& dhd);
    
  public:
    /*! \brief Obtain details from mapped dihedral.
     *  \param mapped the dihedral matched. */
    void MappedWith(Dihedral& mapped);
    
    /*! \brief Apply the parameterisation.
     *  \details Applies the parameteristion. Doing so sets the dihedral type to
     *  the most common type. Additionally, doing so will mean that no more
     *  angles can be mapped with. If the \p self_consistent flag is true, there
     *  is a requirement that all the mapped angles have the same parameters. As
     *  such, if there arises a situtation where this requirement is not met, an
     *  exception will be thrown.
     *  \param self_consistent if the parameterisation needs to be self
     *  consistent. */
    bool ApplyParameterisation(bool self_consistent);
    
  public:
    /*! \brief Calculate the number of mapped dihedrals.
     *  \return the number of dihedrals which have currently been mapped. */
    int64_t NumSourceDihedral() const;
    
    /*! \brief Get the atoms this parameterises the angle between.
     *  \return the three atoms defining the angle this parameterises. */
    stdx::quad<Atom&> GetParameterisedAtoms() const;
    
    /*! \brief Get the dihedral that is parameterised.
     *  \return the parameterised dihedral. */
    Dihedral& GetDihedral() const;
    
    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    TypeGroup GetMostCommonType() const;
    
    const TypeCounts& GetMappedTypeCounts() const;
    
  public:
    bool operator==(const ParamDihedral& dhd) const;
    bool operator!=(const ParamDihedral& dhd) const;
    bool operator<(const  ParamDihedral& dhd) const;
    bool operator>(const  ParamDihedral& dhd) const;
    bool operator<=(const ParamDihedral& dhd) const;
    bool operator>=(const ParamDihedral& dhd) const;
    operator bool() const;
    
  private:
    struct ParamDihedralImpl;
    std::shared_ptr<ParamDihedralImpl> m_pdhddat;
  };
  std::ostream& operator<<(std::ostream& os, const ParamDihedral& dhd);
  
  class ParamMolecule {
  public:
    //! \brief Container type to hold parameterised atoms
    using ParamAtoms = eastl::vector_map<sAtom, ParamAtom>;
    //! \brief Type defining bonds
    using PBond = std::pair<sAtom, sAtom>;
    //! \brief Container type to hold parameterised bonds
    using ParamBonds = eastl::vector_map<PBond, ParamBond>;
    //! \brief Type defining angles
    using PAngle = stdx::triple<sAtom>;
    //! \brief Container type to hold parameterised angles
    using ParamAngles = eastl::vector_map<PAngle, ParamAngle>;
    //! \brief Type defining dihedrals
    using PDihedral = stdx::quad<sAtom>;
    //! \brief Container type to hold parameterised dihedrals
    using ParamDihedrals = eastl::vector_map<PDihedral, ParamDihedral>;
  public:
    ParamMolecule();
    ParamMolecule(const ParamMolecule& mol);
    ParamMolecule(ParamMolecule&& mol);
    ParamMolecule& operator=(const ParamMolecule& mol);
    ParamMolecule& operator=(ParamMolecule&& mol);
    /*! \brief Normal constructor
     *  \param mol the molecule to parameterise. */
    ParamMolecule(Molecule& mol);
    
  public:
    /*! \brief Applies the current parameterisation state.
     *  \details Goes through all parameterised parts and calls the apply
     *  parameterisation method on each of them.
     *  \param self_consistent apply parapmeterisation self consistently. */
    void ApplyParameteristion(bool self_consistent);
    
  public:
    /*! \brief Get the parameterisation of an atom
     *  \param atm the atom to get
     *  \return the parameterisation atom. */
    ParamAtom GetAtom(Atom& atm) const;
    
    /*! \brief Get the parameterisation of a bond
     *  \param bnd the bond to get
     *  \return the parameterisation bond. */
    ParamBond GetBond(Bond& bnd) const;
    
    /*! \brief Get the parameterisation of a bond
     *  \param atms the pair of atoms the bond is between
     *  \return the parameterisation bond. */
    ParamBond GetBond(PBond atms) const;
    
    /*! \brief Get the parameterisation of an angle.
     *  \param ang the angle to get.
     *  \return the parameterisation angle. */
    ParamAngle GetAngle(Angle& ang) const;
    
    /*! \brief Get the parameterisation of an angle.
     *  \param atms the triple of atoms the angle is between.
     *  \return the parameterisation angle. */
    ParamAngle GetAngle(PAngle atms) const;
    
    /*! \brief Get the parameterisation of a dihedral.
     *  \param dhd the dihedral to get.
     *  \return the parameterisation dihedral. */
    ParamDihedral GetDihedral(Dihedral& dhd);
    
    /*! \brief Get the parameterisation of a dihedral.
     *  \details To allow for the parameterisation of dihedrals, such as
     *  impropers used to keep chirality, which are not percieved within the
     *  molecule, if the given atoms are not found, a new ParamDihedral will be
     *  created.
     *  \param atms the quad of atoms the dihedral is between.
     *  \return the parameterisation dihedral. */
    ParamDihedral GetDihedral(PDihedral atms);
    
  private:
    struct ParamMoleculeImpl;
    std::shared_ptr<ParamMoleculeImpl> m_pmoldat;
  };
  
}

#endif /* INDIGOX_CLASSES_PARAMETERISED_HPP */
