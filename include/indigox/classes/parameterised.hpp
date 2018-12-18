#include "../utils/fwd_declares.hpp"
#include "../utils/numerics.hpp"
#include "../utils/quad.hpp"
#include "../utils/triple.hpp"

#include <EASTL/vector_map.h>
#include <array>
#include <vector>

#ifndef INDIGOX_CLASSES_PARAMETERISED_HPP
#define INDIGOX_CLASSES_PARAMETERISED_HPP

namespace indigox {

  class ParamAtom {
  public:
    //! \brief Type giving counts of each mapped atom type
    using TypeCounts = eastl::vector_map<FFAtom, size_t>;
    //! \brief Type giving all mapped charges
    using MappedCharge = std::vector<double>;
    friend class ParamMolecule;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(ParamAtom);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(ParamAtom, atm);

  public:
    /*! \brief Normal constructor
     *  \param atm the Atom this is parameterising. */
    ParamAtom(const Atom &atm);

  public:
    /*! \brief Obtain details from mapped atom
     *  \param mapped the atom matched. */
    void MappedWith(const Atom &mapped);

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
    const Atom &GetAtom() const;

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
    const FFAtom& GetMostCommonType() const;

    const TypeCounts &GetMappedTypeCounts() const;
    const MappedCharge &GetMappedCharges() const;

  private:
    struct ParamAtomImpl;
    std::shared_ptr<ParamAtomImpl> m_data;
  };

  class ParamBond {
  public:
    //! \brief Type giving counts of each mapped bond type
    using TypeCounts = eastl::vector_map<FFBond, size_t>;
    //! \brief Type giving the atoms this parameterised bond is between
    using BondAtoms = std::pair<Atom, Atom>;
    friend class ParamMolecule;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(ParamBond);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(ParamBond, bnd);

  public:
    /*! \brief Normal constructor
     *  \param a,b the atoms which mark the bond this is parameterising. */
    ParamBond(BondAtoms atms, const Bond &bnd);

  public:
    /*! \brief Obtain details from mapped bonds.
     *  \param mapped the bond matched. */
    void MappedWith(const Bond &mapped);

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
    const BondAtoms& GetAtoms() const;

    /*! \brief Get the bond that is parameterised.
     *  \return the parameterised bond. */
    const Bond &GetBond() const;

    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    const FFBond& GetMostCommonType() const;

    const TypeCounts &GetMappedTypeCounts() const;

  private:
    struct ParamBondImpl;
    std::shared_ptr<ParamBondImpl> m_data;
  };

  class ParamAngle {
  public:
    //! \brief Type giving counts of each mapped angle type
    using TypeCounts = eastl::vector_map<FFAngle, size_t>;
    //! \brief Type giving the atoms this parameterised angle is between
    using AngleAtoms = stdx::triple<Atom>;
    friend class ParamMolecule;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(ParamAngle);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(ParamAngle, ang);

  public:
    /*! \brief Normal constructor.
     *  \param a,b,c the atoms which mark the parameterised angle. */
    ParamAngle(AngleAtoms atms, const Angle &ang);

  public:
    /*! \brief Obtain details from mapped angles.
     *  \param mapped the angle matched. */
    void MappedWith(const Angle &mapped);

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
    const AngleAtoms& GetAtoms() const;

    /*! \brief Get the angle that is parameterised.
     *  \return the parameterised bond. */
    const Angle &GetAngle() const;

    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    const FFAngle& GetMostCommonType() const;

    const TypeCounts &GetMappedTypeCounts() const;

  private:
    struct ParamAngleImpl;
    std::shared_ptr<ParamAngleImpl> m_data;
  };

  class ParamDihedral {
  public:
    //! \brief Type giving counts of each mapped dihedral type
    using TypeGroup = std::vector<FFDihedral>;
    using TypeCounts = eastl::vector_map<TypeGroup, size_t>;
    //! \brief Type giving the atoms this parameterised dihedral is between
    using DihedralAtoms = stdx::quad<Atom>;

    friend class ParamMolecule;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(ParamDihedral);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(ParamDihedral, dhd);

  public:
    /*! \brief Normal constructor.
     *  \param a,b,c,d the atoms which mark the parameterised dihedral. */
    ParamDihedral(DihedralAtoms atms, const Dihedral &dhd);

  public:
    /*! \brief Obtain details from mapped dihedral.
     *  \param mapped the dihedral matched. */
    void MappedWith(const Dihedral &mapped);

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
    const DihedralAtoms& GetAtoms() const;

    /*! \brief Get the dihedral that is parameterised.
     *  \return the parameterised dihedral. */
    const Dihedral &GetDihedral() const;

    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    const TypeGroup& GetMostCommonType() const;

    const TypeCounts &GetMappedTypeCounts() const;

  private:
    struct ParamDihedralImpl;
    std::shared_ptr<ParamDihedralImpl> m_data;
  };

  class ParamMolecule {
  public:
    //! \brief Container type to hold parameterised atoms
    using ParamAtoms = eastl::vector_map<Atom, uint32_t>;
    //! \brief Type defining bonds
    using PBond = std::pair<Atom, Atom>;
    //! \brief Container type to hold parameterised bonds
    using ParamBonds = eastl::vector_map<PBond, uint32_t>;
    //! \brief Type defining angles
    using PAngle = stdx::triple<Atom>;
    //! \brief Container type to hold parameterised angles
    using ParamAngles = eastl::vector_map<PAngle, uint32_t>;
    //! \brief Type defining dihedrals
    using PDihedral = stdx::quad<Atom>;
    //! \brief Container type to hold parameterised dihedrals
    using ParamDihedrals = eastl::vector_map<PDihedral, uint32_t>;

  public:
    INDIGOX_GENERIC_PIMPL_CLASS_DEFAULTS(ParamMolecule);
    INDIGOX_GENERIC_PIMPL_CLASS_OPERATORS(ParamMolecule, mol);
    
    /*! \brief Normal constructor
     *  \param mol the molecule to parameterise. */
    ParamMolecule(const Molecule &mol);

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
    const ParamAtom& GetAtom(const Atom &atm) const;

    /*! \brief Get the parameterisation of a bond
     *  \param bnd the bond to get
     *  \return the parameterisation bond. */
    const ParamBond& GetBond(const Bond &bnd) const;

    /*! \brief Get the parameterisation of a bond
     *  \param atms the pair of atoms the bond is between
     *  \return the parameterisation bond. */
    const ParamBond& GetBond(const Atom& a, const Atom& b) const;

    /*! \brief Get the parameterisation of an angle.
     *  \param ang the angle to get.
     *  \return the parameterisation angle. */
    const ParamAngle& GetAngle(const Angle &ang) const;

    /*! \brief Get the parameterisation of an angle.
     *  \param atms the triple of atoms the angle is between.
     *  \return the parameterisation angle. */
    const ParamAngle& GetAngle(const Atom& a, const Atom& b, const Atom& c) const;

    /*! \brief Get the parameterisation of a dihedral.
     *  \param dhd the dihedral to get.
     *  \return the parameterisation dihedral. */
    const ParamDihedral& GetDihedral(const Dihedral &dhd);

    /*! \brief Get the parameterisation of a dihedral.
     *  \details To allow for the parameterisation of dihedrals, such as
     *  impropers used to keep chirality, which are not percieved within the
     *  molecule, if the given atoms are not found, a new ParamDihedral will be
     *  created.
     *  \param atms the quad of atoms the dihedral is between.
     *  \return the parameterisation dihedral. */
    const ParamDihedral& GetDihedral(const Atom& a, const Atom& b, const Atom& c, const Atom& d);

    const std::vector<ParamAtom>& GetAtoms() const;
    const std::vector<ParamBond>& GetBonds() const;
    const std::vector<ParamAngle>& GetAngles() const;
    const std::vector<ParamDihedral>& GetDihedrals() const;

  private:
    struct ParamMoleculeImpl;
    std::shared_ptr<ParamMoleculeImpl> m_data;
  };

} // namespace indigox

#endif /* INDIGOX_CLASSES_PARAMETERISED_HPP */
