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
  
  class IXParamAtom {
  public:
    //! \brief Type giving counts of each mapped atom type
    using TypeCounts = eastl::vector_map<FFAtom, size_>;
    //! \brief Type giving all mapped charges
    using MappedCharge = std::vector<float_>;
    //! \brief Type giving all the mapped atoms
    using MappedAtoms = std::vector<_Atom>;
    
    IXParamAtom() = delete; // no default constructor
    /*! \brief Normal constructor
     *  \param atm the Atom this is parameterising. */
    IXParamAtom(Atom atm) : _atm(atm), _applied(false) { }
    
    /*! \brief Calculates the number of mapped atoms
     *  \return the number of atoms which have currently been mapped to me. */
    size_ NumSourceAtoms() const { return _atms.size(); }
    
    /*! \brief Obtain details from mapped atom
     *  \param mapped the atom matched. */
    void MappedWith(Atom mapped);
    
    /*! \brief Get the atom this parameterises.
     *  \return the atom this parameterises. */
    Atom GetParameterisedAtom() const { return _atm.lock(); }
    
    /*! \brief Get the mean of the charges
     *  \return the mean of the charges. */
    float_ MeanCharge() const {
      return CalculateMean(_charges.begin(), _charges.end()); }
    
    /*! \brief Get the median of the cahrges.
     *  \return the median of the charges. */
    float_ MeadianCharge() {
      return CalculateMedian(_charges.begin(), _charges.end()); }
    
    /*! \brief Get the standard deviation of the charges.
     *  \return the standard deviation of the charges. */
    float_ StandardDeviationCharge() const {
      return CalculateStandardDeviation(_charges.begin(), _charges.end()); }
    
    /*! \brief Get the mode of mapped atom types.
     *  \return the most commonly mapped atom type. */
    FFAtom GetMostCommonType() const {
      if (_counts.empty()) return FFAtom();
      return std::max_element(_counts.begin(), _counts.end(),
                              [](const TypeCounts::value_type& a,
                                 const TypeCounts::value_type& b) {
                                return a.second < b.second;
                              })->first;
    }
    
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
    void ApplyParameterisation(bool self_consistent);
     
    
  private:
    //! \brief Atom which this parameterises
    _Atom _atm;
    //! \brief Mapped atom types
    TypeCounts _counts;
    //! \brief Mapped charges
    MappedCharge _charges;
    //! \brief Mapped atoms
    MappedAtoms _atms;
    //! \brief If parameterisation has already been applied
    bool _applied;
  };
  
  class IXParamBond {
  public:
    //! \brief Type giving counts of each mapped bond type
    using TypeCounts = eastl::vector_map<FFBond, size_>;
    //! \brief Type giving all the mapped bonds
    using MappedBonds = std::vector<_Bond>;
    //! \brief Type giving the atoms this parameterised bond is between
    using BondAtoms = std::pair<_Atom, _Atom>;
    
    IXParamBond() = delete;
    /*! \brief Normal constructor
     *  \param a,b the atoms which mark the bond this is parameterising. */
    IXParamBond(std::pair<Atom, Atom> atms, Bond bnd)
    : _atms(atms), _bnd(bnd), _applied(false) { }
    
    /*! \brief Calculates the number of mapped bonds.
     *  \return the number of bonds which have currently been mapped. */
    size_ NumSourceBonds() const { return _bnds.size(); }
    
    /*! \brief Obtain details from mapped bonds.
     *  \param mapped the bond matched. */
    void MappedWith(Bond mapped);
    
    /*! \brief Get the atoms this parameterises the bond between.
     *  \return the two atoms defining the bond this parameterises. */
    std::pair<Atom, Atom> GetParameterisedAtoms() const {
      return std::make_pair(_atms.first.lock(), _atms.second.lock());
    }
    
    /*! \brief Get the bond that is parameterised.
     *  \return the parameterised bond. */
    Bond GetParameterisedBond() const { return _bnd.lock(); }
    
    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    FFBond GetMostCommonType() const {
      if (_counts.empty()) return FFBond();
      return std::max_element(_counts.begin(), _counts.end(),
                              [](const TypeCounts::value_type& a,
                                 const TypeCounts::value_type& b) {
                                return a.second < b.second;
                              })->first;
    }
    
    /*! \brief Apply the parameterisation.
     *  \details Applies the parameteristion. Doing so sets the bond type to the
     *  most common type. Additionally, doing so will mean that no more bonds
     *  can be mapped with. If the \p self_consistent flag is true, there is a
     *  requirement that all the mapped bonds have the same parameters. As such,
     *  if there arises a situtation where this requirement is not met, an
     *  exception will be thrown.
     *  \param self_consistent if the parameterisation needs to be self
     *  consistent. */
    void ApplyParameterisation(bool self_consistent);
    
  private:
    //! \brief Atoms between which the bond is parameterised
    BondAtoms _atms;
    //! \brief Bond which is parameterised
    _Bond _bnd;
    //! \brief Mapped bond types
    TypeCounts _counts;
    //! \brief Mapped bonds
    MappedBonds _bnds;
    //! \brief If parameterisation has already been applied
    bool _applied;
  };
  
  class IXParamAngle {
  public:
    //! \brief Type giving counts of each mapped angle type
    using TypeCounts = eastl::vector_map<FFAngle, size_>;
    //! \brief Type giving all the mapped angles
    using MappedAngles = std::vector<_Angle>;
    //! \brief Type giving the atoms this parameterised angle is between
    using AngleAtoms = stdx::triple<_Atom, _Atom, _Atom>;
    
    IXParamAngle() = delete;
    /*! \brief Normal constructor.
     *  \param a,b,c the atoms which mark the parameterised angle. */
    IXParamAngle(stdx::triple<Atom, Atom, Atom> atms, Angle ang)
    : _atms(atms), _ang(ang), _applied(false) { }
    
    
    /*! \brief Calculate the number of mapped angles.
     *  \return the number of angles which have currently been mapped. */
    size_ NumSourceAngles() const { return _angs.size(); }
    
    /*! \brief Obtain details from mapped angles.
     *  \param mapped the angle matched. */
    void MappedWith(Angle mapped);
    
    /*! \brief Get the atoms this parameterises the angle between.
     *  \return the three atoms defining the angle this parameterises. */
    stdx::triple<Atom, Atom, Atom> GetParameterisedAtoms() const {
      return stdx::make_triple(_atms.first.lock(),
                               _atms.second.lock(),
                               _atms.third.lock());
    }
    
    /*! \brief Get the angle that is parameterised.
     *  \return the parameterised bond. */
    Angle GetParameterisedAngle() const { return _ang.lock(); }
    
    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    FFAngle GetMostCommonType() const {
      if (_counts.empty()) return FFAngle();
      return std::max_element(_counts.begin(), _counts.end(),
                              [](const TypeCounts::value_type& a,
                                 const TypeCounts::value_type& b) {
                                return a.second < b.second;
                              })->first;
    }
    
    /*! \brief Apply the parameterisation.
     *  \details Applies the parameteristion. Doing so sets the angle type to
     *  the most common type. Additionally, doing so will mean that no more
     *  angles can be mapped with. If the \p self_consistent flag is true, there
     *  is a requirement that all the mapped angles have the same parameters. As
     *  such, if there arises a situtation where this requirement is not met, an
     *  exception will be thrown.
     *  \param self_consistent if the parameterisation needs to be self
     *  consistent. */
    void ApplyParameterisation(bool self_consistent);
    
  private:
    //! \brief Atoms between which the angle is parameterised
    AngleAtoms _atms;
    //! \brief Angle which is parameterised
    _Angle _ang;
    //! \brief Mapped angle types
    TypeCounts _counts;
    //! \brief Mapped angles
    MappedAngles _angs;
    //! \brief If parameterisation has already been applied
    bool _applied;
  };
  
  class IXParamDihedral {
  public:
    //! \brief Type giving counts of each mapped dihedral type
    using TypeCounts = eastl::vector_map<FFDihedral, size_>;
    //! \brief Type giving all the mapped dihedrals
    using MappedDihedrals = std::vector<_Dihedral>;
    //! \brief Type giving the atoms this parameterised dihedral is between
    using DihedralAtoms = stdx::quad<_Atom, _Atom, _Atom, _Atom>;
    
    IXParamDihedral() = delete;
    /*! \brief Normal constructor.
     *  \param a,b,c,d the atoms which mark the parameterised dihedral. */
    IXParamDihedral(stdx::quad<Atom, Atom, Atom, Atom> atms, Dihedral dhd)
    : _atms(atms), _dhd(dhd), _applied(false) { }
    
    /*! \brief Calculate the number of mapped dihedrals.
     *  \return the number of dihedrals which have currently been mapped. */
    size_ NumSourceDihedral() const { return _dhds.size(); }
    
    /*! \brief Obtain details from mapped dihedral.
     *  \param mapped the dihedral matched. */
    void MappedWith(Dihedral mapped);
    
    /*! \brief Get the atoms this parameterises the angle between.
     *  \return the three atoms defining the angle this parameterises. */
    stdx::quad<Atom, Atom, Atom, Atom> GetParameterisedAtoms() const {
      return stdx::make_quad(_atms.first.lock(),
                             _atms.second.lock(),
                             _atms.third.lock(),
                             _atms.fourth.lock());
    }
    
    /*! \brief Get the dihedral that is parameterised.
     *  \return the parameterised dihedral. */
    Dihedral GetParameterisedDihedral() const { return _dhd.lock(); }
    
    /*! \brief Get the mode of mapped bond types.
     *  \return the most commonly mapped bond type. */
    FFDihedral GetMostCommonType() const {
      if (_counts.empty()) return FFDihedral();
      return std::max_element(_counts.begin(), _counts.end(),
                              [](const TypeCounts::value_type& a,
                                 const TypeCounts::value_type& b) {
                                return a.second < b.second;
                              })->first;
    }
    
    /*! \brief Apply the parameterisation.
     *  \details Applies the parameteristion. Doing so sets the dihedral type to
     *  the most common type. Additionally, doing so will mean that no more
     *  angles can be mapped with. If the \p self_consistent flag is true, there
     *  is a requirement that all the mapped angles have the same parameters. As
     *  such, if there arises a situtation where this requirement is not met, an
     *  exception will be thrown.
     *  \param self_consistent if the parameterisation needs to be self
     *  consistent. */
    void ApplyParameterisation(bool self_consistent);
    
  private:
    //! \brief Atoms between which the dihedral is parameterised
    DihedralAtoms _atms;
    //! \brief Dihedral which is parameterised
    _Dihedral _dhd;
    //! \brief Mapped dihedral types
    TypeCounts _counts;
    //! \brief Mapped dihedrals
    MappedDihedrals _dhds;
    //! \brief If parameterisation has already been applied
    bool _applied;
  };
  
  class IXParamMolecule {
  public:
    //! \brief Container type to hold parameterised atoms
    using ParamAtoms = eastl::vector_map<Atom, ParamAtom>;
    //! \brief Type defining bonds
    using PBond = std::pair<Atom, Atom>;
    //! \brief Container type to hold parameterised bonds
    using ParamBonds = eastl::vector_map<PBond, ParamBond>;
    //! \brief Type defining angles
    using PAngle = stdx::triple<Atom, Atom, Atom>;
    //! \brief Container type to hold parameterised angles
    using ParamAngles = eastl::vector_map<PAngle, ParamAngle>;
    //! \brief Type defining dihedrals
    using PDihedral = stdx::quad<Atom, Atom, Atom, Atom>;
    //! \brief Container type to hold parameterised dihedrals
    using ParamDihedrals = eastl::vector_map<PDihedral, ParamDihedral>;
    
    IXParamMolecule() = delete;
    /*! \brief Normal constructor
     *  \param mol the molecule to parameterise. */
    IXParamMolecule(Molecule mol);
    
    /*! \brief Get the parameterisation of an atom
     *  \param atm the atom to get
     *  \return the parameterisation atom. */
    ParamAtom GetAtom(Atom atm) const;
    
    /*! \brief Get the parameterisation of a bond
     *  \param bnd the bond to get
     *  \return the parameterisation bond. */
    ParamBond GetBond(Bond bnd) const;
    
    /*! \brief Get the parameterisation of a bond
     *  \param atms the pair of atoms the bond is between
     *  \return the parameterisation bond. */
    ParamBond GetBond(PBond atms) const;
    
    /*! \brief Get the parameterisation of an angle.
     *  \param ang the angle to get.
     *  \return the parameterisation angle. */
    ParamAngle GetAngle(Angle ang) const;
    
    /*! \brief Get the parameterisation of an angle.
     *  \param atms the triple of atoms the angle is between.
     *  \return the parameterisation angle. */
    ParamAngle GetAngle(PAngle atms) const;
    
    /*! \brief Get the parameterisation of a dihedral.
     *  \param dhd the dihedral to get.
     *  \return the parameterisation dihedral. */
    ParamDihedral GetDihedral(Dihedral dhd);
    
    /*! \brief Get the parameterisation of a dihedral.
     *  \details To allow for the parameterisation of dihedrals, such as
     *  impropers used to keep chirality, which are not percieved within the
     *  molecule, if the given atoms are not found, a new ParamDihedral will be
     *  created.
     *  \param atms the quad of atoms the dihedral is between.
     *  \return the parameterisation dihedral. */
    ParamDihedral GetDihedral(PDihedral atms);
    
    /*! \brief Applies the current parameterisation state.
     *  \details Goes through all parameterised parts and calls the apply
     *  parameterisation method on each of them.
     *  \param self_consistent apply parapmeterisation self consistently. */
    void ApplyParameteristion(bool self_consistent);
    
  private:
    //! \brief Molecule I parameterise
    Molecule _mol;
    //! \brief Parameterised atoms
    ParamAtoms _atms;
    //! \brief Parameterised bonds
    ParamBonds _bnds;
    //! \brief Parameterised angles
    ParamAngles _angs;
    //! \brief Parameterised dihedrals
    ParamDihedrals _dhds;
  };
  
}

#endif /* INDIGOX_CLASSES_PARAMETERISED_HPP */
