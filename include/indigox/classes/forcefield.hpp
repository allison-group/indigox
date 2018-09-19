/*! \file forcefield.hpp */
#include <array>
#include <initializer_list>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

#include <EASTL/bitset.h>

#include "../utils/common.hpp"
#include "../utils/numerics.hpp"
#include "../utils/fwd_declares.hpp"

#ifndef INDIGOX_CLASSES_FORCEFIELD_HPP
#define INDIGOX_CLASSES_FORCEFIELD_HPP

namespace indigox {
  
  //! \brief type for collection of parameter inputs
  using FFParam = std::initializer_list<float_>;
  
  
  class IXFFDihedral {
  public:
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
    //! \brief Friendship allows IXFFDihedral internals to be tested.
    friend struct indigox::test::TestFFDihedral;
    //! \brief Friendship allows dihedral types to be added to a forcefield.
    friend class IXForcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
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
    using DataStore = std::array<float_, STORE_SIZE>;
    //! \brief Type used for storing allowed parameters
    using AllowedMask = eastl::bitset<num_allow_positions>;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFDihedral>& construct,
                                   const uint32_t version);
    
  public:
    IXFFDihedral() = delete; // no default constructor
    
    /*! \brief Normal construtor
     *  \details Creates a new dihedral term of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function.
     *  \param ff the forcefield this type will beloing to. */
    IXFFDihedral(Type type, int_ id, FFParam parameters, const Forcefield& ff);
    
    /*! \brief Empty constructor
     *  \details Creates a new dihedral term of the specified type but without
     *  any parameters set. Used for serialisation purposes.
     *  \param type the type of the potential energy function to use.
     *  \param ff the forcefield to assign to. */
    IXFFDihedral(Type type, const Forcefield& ff);
    
    /*! \brief Get the phase shift (in degrees).
     *  \return the phase shift.
     *  \throws std::runtime_error if the type does not allow a phase shift. */
    inline float_ GetPhaseShift() const {
      return _mask[Allow_PhaseShift] ? _dat[Store_PhaseShift]
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the force constant.
     *  \return the force constant.
     *  \throws std::runtime_error if the type does not allow a force constant.
     */
    inline float_ GetForceConstant() const {
      return _mask[Allow_ForceConstant] ? _dat[Store_ForceConstant]
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the multiplicity.
     *  \return the multiplicity.
     *  \throws std::runtime_error if the type does not allow a multiplicity. */
    inline uint_ GetMultiplicity() const {
      return _mask[Allow_Multiplicity] ? static_cast<uint_>(_dat[Store_Multiplicity])
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the ideal angle (in degrees).
     *  \return the ideal angle.
     *  \throws std::runtime_error if the type does not allow an ideal angle. */
    inline float_ GetIdealAngle() const {
      return _mask[Allow_IdealAngle] ? _dat[Store_IdealAngle]
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the potential energy function type.
     *  \return the type of the potential energy function. */
    inline Type GetType() const { return _type; }
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    inline int_ GetID() const { return _id; }
    
    /*! \brief Get the forcefield this type belongs to.
     *  \return the forcefield. */
    inline Forcefield GetForcefield() const { return _ff.lock(); }
    
  private:
    //! The type of the potential energy function.
    Type _type;
    //! An id for this dihedral type.
    int_ _id;
    //! Terms for the potential energy function.
    DataStore _dat;
    //! Mask for the allowed parameters.
    AllowedMask _mask;
    //! Forcefield
    _Forcefield _ff;
  };
  using DihedralType = IXFFDihedral::Type;
  
  class IXFFAngle {
  public:
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
    //! \brief Friendship allows IXFFAngle internals to be tested.
    friend struct indigox::test::TestFFAngle;
    //! \brief Friendship allows angle types to be added to a forcefield.
    friend class IXForcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
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
    using DataStore = std::array<float_, STORE_SIZE>;
    //! \brief Type used for storing allowed parameters
    using AllowedMask = eastl::bitset<num_allow_positions>;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFAngle>& construct,
                                   const uint32_t version);
    
  public:
    IXFFAngle() = delete;   // no default constructor
    
    /*! \brief Normal constructor
     *  \details Creates a new angle of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function.
     *  \param ff the forcefield this type will belong to. */
    IXFFAngle(Type type, int_ id, FFParam parameters, const Forcefield& ff);
    
    /*! \brief Empty constructor
     *  \details Creates a new dihedral term of the specified type but without
     *  any parameters set. Used for serialisation purposes.
     *  \param type the type of the potential energy function to use.
     *  \param ff the forcefield to assign to. */
    IXFFAngle(Type type, const Forcefield& ff);
    
    /*! \brief Get the force constant.
     *  \return the force constant.
     *  \throws std::runtime_error if the type does not allow a force constant.
     */
    inline float_ GetForceConstant() const {
      return _mask[Allow_ForceConstant] ? _dat[Store_ForceConstant]
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the ideal angle (in degrees).
     *  \return the ideal angle.
     *  \throws std::runtime_error if the type does not allow an ideal angle. */
    inline float_ GetIdealAngle() const {
      return _mask[Allow_IdealAngle] ? _dat[Store_IdealAngle]
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the potential energy function type.
     *  \return the type of the potential energy function. */
    inline Type GetType() const { return _type; }
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    inline int_ GetID() const { return _id; }
    
    /*! \brief Get the linked type for this type.
     *  \details A linked type is an equivalent type with a different functional
     *  type.
     *  \return the linked type. */
    inline FFAngle GetLinkedType() const { return _link; }
    
    /*! \brief Get the forcefield this type belongs to.
     *  \return the forcefield. */
    inline Forcefield GetForcefield() const { return _ff.lock(); }
    
  private:
    //! The type of the potential energy function.
    Type _type;
    //! An id for this angle type.
    int_ _id;
    //! Terms for the potential energy function.
    DataStore _dat;
    //! Mask for the allowed parameters.
    AllowedMask _mask;
    //! Linked angle type
    FFAngle _link;
    //! Forcefield
    _Forcefield _ff;
  };
  using AngleType = IXFFAngle::Type;
  
  class IXFFBond {
  public:
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
    //! \brief Friendship allows IXFFBond internals to be tested.
    friend struct indigox::test::TestFFBond;
    //! \brief Friendship allows bond types to be added to a forcefield.
    friend class IXForcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
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
    
    //! \brief Type used to store all the internal data
    using DataStore = std::array<float_, STORE_SIZE>;
    //! \brief Type used for storing allowed parameters
    using AllowedMask = eastl::bitset<num_allow_positions>;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFBond>& construct,
                                   const uint32_t version);
    
  public:
    IXFFBond() = delete;   // no default constructor
    
    /*! \brief Normal constructor
     *  \details Creates a new bond of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function.
     *  \param ff the forcefield this type will belong to. */
    IXFFBond(Type type, int_ id, FFParam parameters, const Forcefield& ff);
    
    /*! \brief Empty constructor
     *  \details Creates a new bond term of the specified type but without any
     *  parameters set. Used for serialisation purposes.
     *  \param type the ytpe of the potential energy function to use.
     *  \param ff the forcefield to assign to. */
    IXFFBond(Type type, const Forcefield& ff);
    
    /*! \brief Get the force constant.
     *  \return the force constant.
     *  \throws std::runtime_error if the type does not allow a force constant.
     */
    inline float_ GetForceConstant() const {
      return _mask[Allow_ForceConstant] ? _dat[Store_ForceConstant]
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the ideal angle (in degrees).
     *  \return the ideal angle.
     *  \throws std::runtime_error if the type does not allow an ideal angle. */
    inline float_ GetIdealLength() const {
      return _mask[Allow_IdealLength] ? _dat[Store_IdealLength]
      : throw std::runtime_error("Disallowed parameter type requested");
    }
    
    /*! \brief Get the potential energy function type.
     *  \return the type of the potential energy function. */
    inline Type GetType() const { return _type; }
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    inline int_ GetID() const { return _id; }
    
    /*! \brief Get the linked type for this type.
     *  \details A linked type is an equivalent type with a different functional
     *  type.
     *  \return the linked type. */
    inline FFBond GetLinkedType() const { return _link; }
    
    /*! \brief Get the forcefield this type belongs to.
     *  \return the forcefield. */
    inline Forcefield GetForcefield() const { return _ff.lock(); }
    
  private:
    //! The type of the potential energy function.
    Type _type;
    //! An id for this angle type.
    int_ _id;
    //! Terms for the potential energy function.
    DataStore _dat;
    //! Mask for the allowed parameters.
    AllowedMask _mask;
    //! Linked type
    FFBond _link;
    //! Forcefield
    _Forcefield _ff;
  };
  using BondType = IXFFBond::Type;
  
  class IXFFAtom {
    //! \brief Friendship allows IXFFAtomType internals to be tested.
    friend struct indigox::test::TestFFAtom;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    //! \brief Friendship allows atom types to be added to a forcefield.
    friend class IXForcefield;
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFAtom>& construct,
                                   const uint32_t version);
    
  public:
    IXFFAtom() = delete;  // no default constructor
    
    /*! \brief Normal constructor
     *  \details Creates a new atom term.
     *  \param id the ID of the type.
     *  \param name the name of the atom type.
     *  \param the forcefield this type will be a part of. */
    IXFFAtom(int_ id, string_ name, const Forcefield& ff);
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this atom type. */
    inline int_ GetID() const { return _id; }
    
    /*! \brief Get the name for this type
     *  \return the name of this atom type. */
    inline string_ GetName() const { return _name; }
    
    /*! \brief Get the forcefield this type belongs to.
     *  \return the forcefield. */
    inline Forcefield GetForcefield() const { return _ff.lock(); }
    
  private:
    //! ID of this atom type
    int_ _id;
    //! Name of this atom type
    string_ _name;
    //! Forcefield
    _Forcefield _ff;
  };
  
  class IXForcefield : public std::enable_shared_from_this<IXForcefield> {
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
    //! \brief Friendship allows IXForcefield internals to be tested.
    friend struct indigox::test::TestForcefield;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    
    //! \brief Type used to store the atom types.
    using AtomTypes = std::vector<FFAtom>;
    //! \brief Type used to store the bond types.
    using BondTypes = std::map<BondType, std::vector<FFBond>>;
    //! \brief Type used to store the angle types.
    using AngleTypes = std::map<AngleType, std::vector<FFAngle>>;
    //! \brief Type used to store the dihedral types.
    using DihedralTypes = std::map<DihedralType, std::vector<FFDihedral>>;
    
  public:
    IXForcefield() = delete;  // no default constructor
    
    /*! \brief Normal constructor
     *  \details Create a new forcefield from the given family. The family
     *  dictates which types of functional forms are allowed. Currently
     *  implemented is the GROMOS family of forcefields, which allows harmonic
     *  and quartic bond types, harmonic and caosine-harmonic angle types, and
     *  proper and improper dihedral types.
     *  \param family the family of forcefields to use.
     *  \param name the name of the forcefield. */
    IXForcefield(Family family, string_ name);
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXForcefield>& construct,
                                   const uint32_t version);
    
  private:
    /*! \brief Adds a new bond type to the forcefield.
     *  \param type the type of the bond.
     *  \param id the id of the bond.
     *  \param parameters the parameters of the bond.
     *  \return the newly created bond type. */
    FFBond NewBondType(BondType type, int_ id, FFParam parameters);
    
    /*! \brief Adds a new angle type to the forcefield.
     *  \param type the type of the angle.
     *  \param id the id of the angle.
     *  \param parameters the parameters of the angle.
     *  \return the newly created angle type. */
    FFAngle NewAngleType(AngleType type, int_ id, FFParam parameters);
    
    /*! \brief Adds a new dihedral type to the forcefield.
     *  \param type the type of the dihedral.
     *  \param id the id of the dihedral.
     *  \param parameters the parameters of the dihedral.
     *  \return the newly created dihedral type. */
    FFDihedral NewDihedralType(DihedralType type, int_ id, FFParam parameters);
    
  public:
    //  Atom type handlings
    /*! \brief Add a new atom type to the forcefield
     *  \param id the id of the atom type.
     *  \param name the name of the atom type.
     *  \return the newly created atom type.
     *  \throws std::runtime_error if an atomtype with the given name or id
     *  already exists within the forcefield. */
    FFAtom NewAtomType(int_ id, string_ name);
    
    /*! \brief Get the AtomType with the given name.
     *  \param name the name of the atom type.
     *  \return the atomtype with the name, or an empty shared_ptr. */
    FFAtom GetAtomType(string_ name) const;
    
    /*! \brief Get the AtomType with the given id.
     *  \param id the id of the atom type.
     *  \return the atomtype with the id, or an empty shared_ptr. */
    FFAtom GetAtomType(int_ id) const;
    
    /*! \brief Reserve space for atom types.
     *  \param size the amount of space to reserve. */
    void ReserveAtomTypes(size_ size) { _atms.reserve(size); }
    
    /*! \brief Get the number of atom types in the forcefield.
     *  \return the number of atom types in the forcefield. */
    size_ NumAtomTypes() const { return _atms.size(); }
    
    // Bond type handlings
    /*! \brief Adds a new bond type to the forcefield with two parameters.
     *  \details This constructor should be used for BondType::Harmonic and
     *  BondType::Quartic bond types. In both cases, \p a is the force constant
     *  and \p b is the ideal length.
     *  \param type the type of the bond.
     *  \param id the id of the bond.
     *  \param a,b the two parameter values.
     *  \return the newly created bond type. */
    inline FFBond NewBondType(BondType type, int_ id, float_ a, float_ b) {
      return NewBondType(type, id, {a, b});
    }
    
    /*! \brief Get the type of bond with the given id.
     *  \param type the type of bond to get.
     *  \param the id of the bond type.
     *  \return the requested bond type, or an empty shared_ptr. */
    FFBond GetBondType(BondType type, int_ id) const;
    
    /*! \brief Get the first bond type with the given id.
     *  \details Searches through all allowed bond types in order and returns
     *  the first type with the id.
     *  \param id the id of the bond type.
     *  \return the bond type, or an empty shared_ptr. */
    FFBond GetBondType(int_ id) const;
    
    /*! \brief Link together two equivalent bond types.
     *  \details Two bond types can be linked together when they are regarded as
     *  the same type with different functional forms. For example, the GROMOS
     *  forcefield provides both harmonic and quartic force constants for all
     *  bond types. This method will remove any existing links.
     *  \param a,b the two types to link. */
    void LinkBondTypes(FFBond a, FFBond b);
    
    /*! \brief Reserve space for bond types.
     *  \param type the type of bond to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveBondTypes(BondType type, size_ size) {
      if (_bnds.find(type) != _bnds.end()) _bnds.at(type).reserve(size);
    }
    
    /*! \brief Get the total number of bond types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of bond types in the forcefield. */
    size_ NumBondTypes() const {
      return std::accumulate(_bnds.begin(), _bnds.end(), 0,
                             [](size_ a, auto& b){return a + b.second.size();});
    }
    
    /*! \brief Get the number of bond types of the given functional form.
     *  \details Determines how many of the given bond type funcitonal form are
     *  present in the forcefield. In the case that the provided type is not
     *  supported by the forcefield, 0 is returned.
     *  \param type the type of bond to count.
     *  \return the number of bond types in the forcefield. */
    size_ NumBondTypes(BondType type) const {
      return (_bnds.find(type) != _bnds.end()) ? _bnds.at(type).size() : 0;
    }
    
    // Angle type handlings
    /*! \brief Adds a new angle type to the forcefield with two parameters.
     *  \details This constructor should be used for AngleType::Harmonic and
     *  AngleType::CosineHarmonic bond types. In both cases, \p a is the force
     *  constant and \p b is the ideal angle.
     *  \param type the type of the angle.
     *  \param id the id of the angle.
     *  \param a,b the two parameter values.
     *  \return the newly created angle type. */
    inline FFAngle NewAngleType(AngleType type, int_ id, float_ a, float_ b) {
      return NewAngleType(type, id, {a, b});
    }
    
    /*! \brief Get the type of angle with the given id.
     *  \param type the type of angle to get.
     *  \param the id of the angle type.
     *  \return the requested angle type, or an empty shared_ptr. */
    FFAngle GetAngleType(AngleType type, int_ id) const;
    
    /*! \brief Get the first angle type with the given id.
     *  \details Searches through all allowed angle types in order and returns
     *  the first type with the given id.
     *  \param id the id of the angle type.
     *  \return the angle type, or an empty shared_ptr. */
    FFAngle GetAngleType(int_ id) const;
    
    /*! \brief Link together two equivalent angle types.
     *  \details Two angle types can be linked together when they are regarded
     *  as the same type with different functional forms. For example, the
     *  GROMOS forcefield provides both harmonic and cosine-harmonic force
     *  constants for all angle types. This method will remove any existing
     *  links.
     *  \param a,b the two types to link. */
    void LinkAngleTypes(FFAngle a, FFAngle b);
    
    /*! \brief Reserve space for angle types.
     *  \param type the type of angle to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveAngleTypes(AngleType type, size_ size) {
      if (_angs.find(type) != _angs.end()) _angs.at(type).reserve(size);
    }
    
    /*! \brief Get the number of angle types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of angle types in the forcefield. */
    size_ NumAngleTypes() const {
      return std::accumulate(_angs.begin(), _angs.end(), 0,
                             [](size_ a, auto& b){return a + b.second.size();});
    }
    
    /*! \brief Get the number of angle types of the given functional form.
     *  \details Determines how may of the given angle types are present in the
     *  forcefield. In the case that the provided type is not supposrted by the
     *  forcefield, 0 is returned.
     *  \param type the type of angle to count.
     *  \return the number of angle types in the forcefield. */
    size_ NumAngleTypes(AngleType type) const {
      return (_angs.find(type) != _angs.end()) ? _angs.at(type).size() : 0;
    }
    
    // Dihedral type handlings
    /*! \brief Adds a new dihedral type to the forcefield with three parameters.
     *  \details This constructor should be used for DihedralType::Proper
     *  dihedral types. In that case, a is the force constant, b is the phase
     *  and c is the multiplicity.
     *  \param type the type of dihedral.
     *  \param id the id of the dihedral.
     *  \param a,b,c the three parameter values.
     *  \return the newly created dihedral type. */
    inline FFDihedral NewDihedralType(DihedralType type, int_ id, float_ a,
                                      float_ b, float_ c) {
      return NewDihedralType(type, id, {b,a,c});
    }
    
    /*! \brief Adds a new dihedral type to the forcefield with two parameters.
     *  \details This constructor should be used for DihedralType::Improper
     *  dihedral types. In that case, a is the force constant and b is the ideal
     *  dihedral angle
     *  \param type the type of dihedral.
     *  \param id the id of the dihedral.
     *  \param a,b the three parameter values.
     *  \return the newly created dihedral type. */
    inline FFDihedral NewDihedralType(DihedralType type, int_ id, float_ a,
                                      float_ b) {
      return NewDihedralType(type, id, {b,a});
    }
    
    /*! \brief Get the type of dihedral with the given id.
     *  \param type the type of dihedral to get.
     *  \param the id of the dihedral type.
     *  \return the requested dihedral type, or an empty shared_ptr. */
    FFDihedral GetDihedralType(DihedralType type, int_ id) const;
    
    /*! \brief Get the first dihedral type with the given id.
     *  \details Searches through all allowed dihedral types in order and
     *  returns the first type with the given id.
     *  \param id the id of the dihedral type.
     *  \return the dihedral type, or an empty shared_ptr. */
    FFDihedral GetDihedralType(int_ id) const;
    
    /*! \brief Reserve space for dihedral types.
     *  \param type the type of dihedral to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveDihedralTypes(DihedralType type, size_ size) {
      if (_dhds.find(type) != _dhds.end()) _dhds.at(type).reserve(size);
    }
    
    /*! \brief Get the number of dihedral types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of bond types in the forcefield. */
    size_ NumDihedralTypes() const {
      return std::accumulate(_dhds.begin(), _dhds.end(), 0,
                             [](auto a, auto& b){return a + b.second.size();});
    }
    
    /*! \brief Get the number of dihedral types of the given functional form.
     *  \details Determines how many of the given dihedral types are present in
     *  the forcefield. In the case that the provided type is not supported by
     *  the forcefield, 0 is returned.
     *  \param type the type of dihedral to count.
     *  \return the number of dihedral types in the forcefield. */
    size_ NumDihedralTypes(DihedralType type) const {
      return (_dhds.find(type) != _dhds.end()) ? _dhds.at(type).size() : 0;
    }
    
    // Normal handlings
    /*! \brief Get the family of the forcefield
     *  \return the forcefields family. */
    inline Family GetFamily() const { return _family; }
    
    /*! \brief Get the name of the forcefield
     *  \return the name of the forcefield. */
    inline string_ GetName() const { return _name; }
    
  private:
    //! The family of the forcefield
    Family _family;
    //! The name of the forcefield
    string_ _name;
    //! Atom types
    AtomTypes _atms;
    //! Bond types
    BondTypes _bnds;
    //! Angle types
    AngleTypes _angs;
    //! Dihedral types;
    DihedralTypes _dhds;
  };
  using FFFamily = IXForcefield::Family;
  
  Forcefield GenerateGROMOS54A7();
}

#endif    /* INDIGOX_CLASSES_FORCEFIELD_HPP */
