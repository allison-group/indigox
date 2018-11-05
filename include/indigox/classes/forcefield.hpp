/*! \file forcefield.hpp */
#include <array>
#include <initializer_list>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <EASTL/bitset.h>

#include "periodictable.hpp"
#include "../utils/common.hpp"
#include "../utils/fwd_declares.hpp"

#ifndef INDIGOX_CLASSES_FORCEFIELD_HPP
#define INDIGOX_CLASSES_FORCEFIELD_HPP

namespace indigox {
  
  //! \brief type for collection of parameter inputs
  using FFParam = std::initializer_list<double>;
  
  class FFAtom {
    //! \brief Friendship allows IXFFAtomType internals to be tested.
    friend struct indigox::test::TestFFAtom;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    //! \brief Friendship allows atom types to be added to a forcefield.
    friend class Forcefield;
    
  private:
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    FFAtom();
    FFAtom(const FFAtom& atm);
    FFAtom(FFAtom&& atm);
    FFAtom& operator=(const FFAtom& atm);
    FFAtom& operator=(FFAtom&& atm);
    
  private:
    /*! \brief Normal constructor
     *  \details Creates a new atom term.
     *  \param id the ID of the type.
     *  \param name the name of the atom type.
     *  \param ff the forcefield this type will be a part of. */
    FFAtom(int32_t id, std::string name, const Element& element,
           uint32_t impH, const Forcefield& ff);
    
  public:
    /*! \brief Get the ID for this type.
     *  \return the ID of this atom type. */
    int32_t GetID() const;
    
    /*! \brief Get the name for this type
     *  \return the name of this atom type. */
    std::string GetName() const;
    
    Forcefield& GetForcefield() const;
    uint32_t GetImplicitHydrogenCount() const;
    Element GetElement() const;
    
  public:
    bool operator==(const FFAtom& atm) const;
    bool operator!=(const FFAtom& atm) const;
    bool operator<(const FFAtom& atm) const;
    bool operator>(const FFAtom& atm) const;
    bool operator<=(const FFAtom& atm) const;
    bool operator>=(const FFAtom& atm) const;
    operator bool() const;
    
  private:
    struct FFAtomImpl;
    std::shared_ptr<FFAtomImpl> m_ffatmdat;
  };
  
  class FFBond {
  public:
    //! \brief Friendship allows IXFFBond internals to be tested.
    friend struct indigox::test::TestFFBond;
    //! \brief Friendship allows bond types to be added to a forcefield.
    friend class Forcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
    /*! \brief Enum for the different function types a bond term can have.
     *  \details Though other values are provided, currently only the harmonic and
     *  quartic function types are technically supported. */
    enum class Type {
      Empty,    //!< An empty bond type with no parameters
      Harmonic, //!< A harmonic bond type taking the form \f$\frac{k_b}{2}(x - x_0)^2\f$.
      Quartic,  //!< A quartic bond type taking the form \f$\frac{k_b}{2}(x^2 - {x_0}^2)^2\f$.
      Morse,
      Cubic,
      FENE
    };
    
  private:
    //! \brief Enum giving the allow positions of each parameter type.
    enum AllowEnum {
      Allow_ForceConstant,
      Allow_IdealLength,
      num_allow_positions
    };
    
    //! \brief Enum giving the stored index of each parameter type.
    enum StoreEnum {
      Store_ForceConstant = 0,
      Store_IdealLength = 1,
      STORE_SIZE = 2
    };
    
    //! \brief Type used for storing allowed parameters
    using AllowedMask = eastl::bitset<num_allow_positions, uint8_t>;
    
    //! \brief Type used to store all the internal data
    using DataStore = std::array<double, STORE_SIZE>;
    
  private:
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    FFBond();
    FFBond(const FFBond& bnd);
    FFBond(FFBond&& bnd);
    FFBond& operator=(const FFBond& bnd);
    FFBond& operator=(FFBond&& bnd);
    
  private:
    /*! \brief Normal constructor
     *  \details Creates a new bond of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function.
     *  \param ff the forcefield this type will belong to. */
    FFBond(Type type, int32_t id, FFParam params, const Forcefield& ff);
    
  public:
    /*! \brief Get the force constant.
     *  \return the force constant.
     *  \throws std::runtime_error if the type does not allow a force constant.
     */
    double GetForceConstant() const;
    
    /*! \brief Get the ideal angle (in degrees).
     *  \return the ideal angle.
     *  \throws std::runtime_error if the type does not allow an ideal angle. */
    double GetIdealLength() const;
    
    /*! \brief Get the potential energy function type.
     *  \return the type of the potential energy function. */
    Type GetType() const;
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    int32_t GetID() const;
    
    /*! \brief Get the linked type for this type.
     *  \details A linked type is an equivalent type with a different functional
     *  type.
     *  \return the linked type. */
    FFBond& GetLinkedType() const;
    
    Forcefield& GetForcefield() const;
    
  public:
    bool operator==(const FFBond& bnd) const;
    bool operator!=(const FFBond& bnd) const;
    bool operator<(const FFBond& bnd) const;
    bool operator>(const FFBond& bnd) const;
    bool operator<=(const FFBond& bnd) const;
    bool operator>=(const FFBond& bnd) const;
    operator bool() const;
    
  private:
    struct FFBondImpl;
    std::shared_ptr<FFBondImpl> m_ffbnddat;
  };
  using BondType = FFBond::Type;
  
  class FFAngle {
  public:
    //! \brief Friendship allows IXFFAngle internals to be tested.
    friend struct indigox::test::TestFFAngle;
    //! \brief Friendship allows angle types to be added to a forcefield.
    friend class Forcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
    /*! \brief Enum for the different function types a angle term can have.
     *  \details Though other values are provided, currently only the harmonic and
     *  cosine-harmonic function types are technically supported. */
    enum class Type {
      Empty,
      Harmonic,
      CosineHarmonic,
      UreyBradley,
      Quartic
    };
    
  private:
    //! \brief Enum giving the allow positions of each parameter type.
    enum AllowEnum {
      Allow_ForceConstant,
      Allow_IdealAngle,
      num_allow_positions
    };
    
    //! \brief Enum giving the stored index of each parameter type.
    enum StoreEnum {
      Store_ForceConstant = 0,
      Store_IdealAngle = 1,
      STORE_SIZE = 2
    };
    
    //! \brief Type used to store all the internal data
    using DataStore = std::array<double, STORE_SIZE>;
    //! \brief Type used for storing allowed parameters
    using AllowedMask = eastl::bitset<num_allow_positions, uint8_t>;
    
  private:
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    FFAngle();
    FFAngle(const FFAngle& ang);
    FFAngle(FFAngle&& ang);
    FFAngle& operator=(const FFAngle& ang);
    FFAngle& operator=(FFAngle&& ang);
    
  private:
    /*! \brief Normal constructor
     *  \details Creates a new angle of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function.
     *  \param ff the forcefield this type will belong to. */
    FFAngle(Type type, int32_t id, FFParam parameters, const Forcefield& ff);
    
  public:
    /*! \brief Get the force constant.
     *  \return the force constant.
     *  \throws std::runtime_error if the type does not allow a force constant.
     */
    double GetForceConstant() const;
    
    /*! \brief Get the ideal angle (in degrees).
     *  \return the ideal angle.
     *  \throws std::runtime_error if the type does not allow an ideal angle. */
    double GetIdealAngle() const;
    
    /*! \brief Get the potential energy function type.
     *  \return the type of the potential energy function. */
    Type GetType() const;
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    int32_t GetID() const;
    
    /*! \brief Get the linked type for this type.
     *  \details A linked type is an equivalent type with a different functional
     *  type.
     *  \return the linked type. */
    FFAngle& GetLinkedType() const;
    
    Forcefield& GetForcefield() const;
    
  public:
    bool operator==(const FFAngle& ang) const;
    bool operator!=(const FFAngle& ang) const;
    bool operator<(const  FFAngle& ang) const;
    bool operator>(const  FFAngle& ang) const;
    bool operator<=(const FFAngle& ang) const;
    bool operator>=(const FFAngle& ang) const;
    operator bool() const;
    
  private:
    struct FFAngleImpl;
    std::shared_ptr<FFAngleImpl> m_ffangdat;
  };
  using AngleType = FFAngle::Type;
  
  class FFDihedral {
  public:
    //! \brief Friendship allows IXFFDihedral internals to be tested.
    friend struct indigox::test::TestFFDihedral;
    //! \brief Friendship allows dihedral types to be added to a forcefield.
    friend class Forcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
    /*! \brief Enum for the different function types a dihedral term can have.
     *  \details Though other values are provided, currently only the proper and
     *  improper function types are technically supported. */
    enum class Type {
      Empty,
      Proper,
      Improper,
      RyckaertBellemans,
      PeriodicImproper,
      Fourier,
      Restricted
    };
    
  private:
    //! \brief Enum giving the allow positions of each parameter type.
    enum AllowEnum {
      Allow_PhaseShift,
      Allow_ForceConstant,
      Allow_Multiplicity,
      Allow_IdealAngle,
      num_allow_positions
    };
    
    //! \brief Enum giving the stored index of each parameter type.
    enum StoreEnum {
      Store_PhaseShift = 0,
      Store_ForceConstant = 1,
      Store_Multiplicity = 2,
      Store_IdealAngle = 0,
      STORE_SIZE = 3
    };
    
    //! \brief Type used to store all the internal data
    using DataStore = std::array<double, STORE_SIZE>;
    //! \brief Type used for storing allowed parameters
    using AllowedMask = eastl::bitset<num_allow_positions>;
    
  private:
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  public:
    FFDihedral();
    FFDihedral(const FFDihedral& dhd);
    FFDihedral(FFDihedral&& dhd);
    FFDihedral& operator=(const FFDihedral& dhd);
    FFDihedral& operator=(FFDihedral&& dhd);
    
  private:
    /*! \brief Normal construtor
     *  \details Creates a new dihedral term of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function.
     *  \param ff the forcefield this type will beloing to. */
    FFDihedral(Type type, int32_t id, FFParam parameters, const Forcefield& ff);
    
  public:
    /*! \brief Get the phase shift (in degrees).
     *  \return the phase shift.
     *  \throws std::runtime_error if the type does not allow a phase shift. */
    double GetPhaseShift() const;
    
    /*! \brief Get the force constant.
     *  \return the force constant.
     *  \throws std::runtime_error if the type does not allow a force constant.
     */
    double GetForceConstant() const;
    
    /*! \brief Get the multiplicity.
     *  \return the multiplicity.
     *  \throws std::runtime_error if the type does not allow a multiplicity. */
    int32_t GetMultiplicity() const;
    
    /*! \brief Get the ideal angle (in degrees).
     *  \return the ideal angle.
     *  \throws std::runtime_error if the type does not allow an ideal angle. */
    double GetIdealAngle() const;
    
    /*! \brief Get the potential energy function type.
     *  \return the type of the potential energy function. */
    Type GetType() const;
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    int32_t GetID() const;
    
    Forcefield& GetForcefield() const;
    
  public:
    bool operator==(const FFDihedral& dhd) const;
    bool operator!=(const FFDihedral& dhd) const;
    bool operator<(const  FFDihedral& dhd) const;
    bool operator>(const  FFDihedral& dhd) const;
    bool operator<=(const FFDihedral& dhd) const;
    bool operator>=(const FFDihedral& dhd) const;
    operator bool() const;
    
  private:
    struct FFDihedralImpl;
    std::shared_ptr<FFDihedralImpl> m_ffdhddat;
  };
  using DihedralType = FFDihedral::Type;
  
  class Forcefield {
    //! \brief Friendship allows IXForcefield internals to be tested.
    friend struct indigox::test::TestForcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
  public:
    /*! \brief Enum for the different families of forcefields.
     *  \details Though other values are provided, currently only the GROMOS
     *  family of forcefields are technically supported. */
    enum class Family {
      Empty,
      GROMOS,
      CHARMM,
      AMBER,
      Other
    };
    
  private:
    //! \brief Type used to store the atom types.
    using AtomTypes = std::vector<FFAtom>;
    //! \brief Type used to store the bond types.
    using BondTypes = std::map<BondType, std::vector<FFBond>>;
    //! \brief Type used to store the angle types.
    using AngleTypes = std::map<AngleType, std::vector<FFAngle>>;
    //! \brief Type used to store the dihedral types.
    using DihedralTypes = std::map<DihedralType, std::vector<FFDihedral>>;
    
  public:
    Forcefield();
    Forcefield(const Forcefield& ff);
    Forcefield(Forcefield&& ff);
    Forcefield& operator=(const Forcefield& ff);
    Forcefield& operator=(Forcefield&& ff);
    
    /*! \brief Normal constructor
     *  \details Create a new forcefield from the given family. The family
     *  dictates which types of functional forms are allowed. Currently
     *  implemented is the GROMOS family of forcefields, which allows harmonic
     *  and quartic bond types, harmonic and caosine-harmonic angle types, and
     *  proper and improper dihedral types.
     *  \param family the family of forcefields to use.
     *  \param name the name of the forcefield. */
    Forcefield(Family family, std::string name);
    
  private:
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t version);
    
  private:
    /*! \brief Adds a new bond type to the forcefield.
     *  \param type the type of the bond.
     *  \param id the id of the bond.
     *  \param parameters the parameters of the bond.
     *  \return the newly created bond type. */
    FFBond& NewBondType(BondType type, int32_t id, FFParam param);
    
    /*! \brief Adds a new angle type to the forcefield.
     *  \param type the type of the angle.
     *  \param id the id of the angle.
     *  \param parameters the parameters of the angle.
     *  \return the newly created angle type. */
    FFAngle& NewAngleType(AngleType type, int32_t id, FFParam param);
    
    /*! \brief Adds a new dihedral type to the forcefield.
     *  \param type the type of the dihedral.
     *  \param id the id of the dihedral.
     *  \param parameters the parameters of the dihedral.
     *  \return the newly created dihedral type. */
    FFDihedral& NewDihedralType(DihedralType type, int32_t id, FFParam param);
    
  public:
    /*! \brief Add a new atom type to the forcefield
     *  \param id the id of the atom type.
     *  \param name the name of the atom type.
     *  \return the newly created atom type.
     *  \throws std::runtime_error if an atomtype with the given name or id
     *  already exists within the forcefield. */
    FFAtom& NewAtomType(int32_t id, std::string name, const Element& element,
                        uint32_t implictH = 0);
    
    /*! \brief Reserve space for atom types.
     *  \param size the amount of space to reserve. */
    void ReserveAtomTypes(size_t size);
    
    /*! \brief Adds a new bond type to the forcefield with two parameters.
     *  \details This constructor should be used for BondType::Harmonic and
     *  BondType::Quartic bond types. In both cases, \p a is the force constant
     *  and \p b is the ideal length.
     *  \param type the type of the bond.
     *  \param id the id of the bond.
     *  \param a,b the two parameter values.
     *  \return the newly created bond type. */
    FFBond& NewBondType(BondType type, int32_t id, double a, double b);
    
    /*! \brief Link together two equivalent bond types.
     *  \details Two bond types can be linked together when they are regarded as
     *  the same type with different functional forms. For example, the GROMOS
     *  forcefield provides both harmonic and quartic force constants for all
     *  bond types. This method will remove any existing links.
     *  \param a,b the two types to link. */
    void LinkBondTypes(FFBond& a, FFBond& b);
    
    /*! \brief Reserve space for bond types.
     *  \param type the type of bond to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveBondTypes(BondType type, size_t size);
    
    /*! \brief Adds a new angle type to the forcefield with two parameters.
     *  \details This constructor should be used for AngleType::Harmonic and
     *  AngleType::CosineHarmonic bond types. In both cases, \p a is the force
     *  constant and \p b is the ideal angle.
     *  \param type the type of the angle.
     *  \param id the id of the angle.
     *  \param a,b the two parameter values.
     *  \return the newly created angle type. */
    FFAngle& NewAngleType(AngleType type, int32_t id, double a, double b);
    
    /*! \brief Link together two equivalent angle types.
     *  \details Two angle types can be linked together when they are regarded
     *  as the same type with different functional forms. For example, the
     *  GROMOS forcefield provides both harmonic and cosine-harmonic force
     *  constants for all angle types. This method will remove any existing
     *  links.
     *  \param a,b the two types to link. */
    void LinkAngleTypes(FFAngle& a, FFAngle& b);
    
    /*! \brief Reserve space for angle types.
     *  \param type the type of angle to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveAngleTypes(AngleType type, size_t size);
    
    /*! \brief Adds a new dihedral type to the forcefield with three parameters.
     *  \details This constructor should be used for DihedralType::Proper
     *  dihedral types. In that case, a is the force constant, b is the phase
     *  and c is the multiplicity.
     *  \param type the type of dihedral.
     *  \param id the id of the dihedral.
     *  \param a,b,c the three parameter values.
     *  \return the newly created dihedral type. */
    FFDihedral& NewDihedralType(DihedralType type, int32_t id, double a,
                                double b, double c);
    
    /*! \brief Adds a new dihedral type to the forcefield with two parameters.
     *  \details This constructor should be used for DihedralType::Improper
     *  dihedral types. In that case, a is the force constant and b is the ideal
     *  dihedral angle
     *  \param type the type of dihedral.
     *  \param id the id of the dihedral.
     *  \param a,b the three parameter values.
     *  \return the newly created dihedral type. */
    FFDihedral& NewDihedralType(DihedralType type, int32_t id, double a,
                                double b);
    
    /*! \brief Reserve space for dihedral types.
     *  \param type the type of dihedral to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveDihedralTypes(DihedralType type, size_t size);
    
  public:
    /*! \brief Get the AtomType with the given name.
     *  \param name the name of the atom type.
     *  \return the atomtype with the name, or an empty shared_ptr. */
    FFAtom& GetAtomType(std::string name) const;
    
    /*! \brief Get the AtomType with the given id.
     *  \param id the id of the atom type.
     *  \return the atomtype with the id, or an empty shared_ptr. */
    FFAtom& GetAtomType(int32_t id) const;
    
    /*! \brief Get the number of atom types in the forcefield.
     *  \return the number of atom types in the forcefield. */
    size_t NumAtomTypes() const;
    
    /*! \brief Get the type of bond with the given id.
     *  \param type the type of bond to get.
     *  \param the id of the bond type.
     *  \return the requested bond type, or an empty shared_ptr. */
    FFBond& GetBondType(BondType type, int32_t id) const;
    
    /*! \brief Get the first bond type with the given id.
     *  \details Searches through all allowed bond types in order and returns
     *  the first type with the id.
     *  \param id the id of the bond type.
     *  \return the bond type, or an empty shared_ptr. */
    FFBond& GetBondType(int32_t id) const;
    
    /*! \brief Get the total number of bond types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of bond types in the forcefield. */
    size_t NumBondTypes() const;
    
    /*! \brief Get the number of bond types of the given functional form.
     *  \details Determines how many of the given bond type funcitonal form are
     *  present in the forcefield. In the case that the provided type is not
     *  supported by the forcefield, 0 is returned.
     *  \param type the type of bond to count.
     *  \return the number of bond types in the forcefield. */
    size_t NumBondTypes(BondType type) const;
    
    /*! \brief Get the type of angle with the given id.
     *  \param type the type of angle to get.
     *  \param the id of the angle type.
     *  \return the requested angle type, or an empty shared_ptr. */
    FFAngle& GetAngleType(AngleType type, int32_t id) const;
    
    /*! \brief Get the first angle type with the given id.
     *  \details Searches through all allowed angle types in order and returns
     *  the first type with the given id.
     *  \param id the id of the angle type.
     *  \return the angle type, or an empty shared_ptr. */
    FFAngle& GetAngleType(int32_t id) const;
    
    /*! \brief Get the number of angle types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of angle types in the forcefield. */
    size_t NumAngleTypes() const;
    
    /*! \brief Get the number of angle types of the given functional form.
     *  \details Determines how may of the given angle types are present in the
     *  forcefield. In the case that the provided type is not supposrted by the
     *  forcefield, 0 is returned.
     *  \param type the type of angle to count.
     *  \return the number of angle types in the forcefield. */
    size_t NumAngleTypes(AngleType type) const;
    
    /*! \brief Get the type of dihedral with the given id.
     *  \param type the type of dihedral to get.
     *  \param the id of the dihedral type.
     *  \return the requested dihedral type, or an empty shared_ptr. */
    FFDihedral& GetDihedralType(DihedralType type, int32_t id) const;
    
    /*! \brief Get the first dihedral type with the given id.
     *  \details Searches through all allowed dihedral types in order and
     *  returns the first type with the given id.
     *  \param id the id of the dihedral type.
     *  \return the dihedral type, or an empty shared_ptr. */
    FFDihedral& GetDihedralType(int32_t id) const;
    
    /*! \brief Get the number of dihedral types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of bond types in the forcefield. */
    size_t NumDihedralTypes() const;
    
    /*! \brief Get the number of dihedral types of the given functional form.
     *  \details Determines how many of the given dihedral types are present in
     *  the forcefield. In the case that the provided type is not supported by
     *  the forcefield, 0 is returned.
     *  \param type the type of dihedral to count.
     *  \return the number of dihedral types in the forcefield. */
    size_t NumDihedralTypes(DihedralType type) const;
    
    /*! \brief Get the family of the forcefield
     *  \return the forcefields family. */
    Family GetFamily() const;
    
    /*! \brief Get the name of the forcefield
     *  \return the name of the forcefield. */
    std::string GetName() const;
    
  public:
    bool operator==(const Forcefield& ff) const;
    bool operator!=(const Forcefield& ff) const;
    operator bool() const;
    
  private:
    struct ForcefieldImpl;
    std::shared_ptr<ForcefieldImpl> m_ffdat;
  };
  using FFFamily = Forcefield::Family;
  
  Forcefield GenerateGROMOS54A7();
}

#endif    /* INDIGOX_CLASSES_FORCEFIELD_HPP */
