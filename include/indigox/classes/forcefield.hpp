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

#ifndef INDIGOX_CLASSES_FORCEFIELD_HPP
#define INDIGOX_CLASSES_FORCEFIELD_HPP

namespace indigox {
  
  class IXForcefield;
  class IXFFAtom;
  class IXFFBond;
  class IXFFAngle;
  class IXFFDihedral;
  namespace test {
    struct TestForcefield;
    struct TestFFAtom;
    struct TestFFBond;
    struct TestFFAngle;
    struct TestFFDihedral;
  }
  
  //! \brief shared_ptr for normal use of the IXForcefield class.
  using Forcefield = std::shared_ptr<IXForcefield>;
  //! \brief shared_ptr for normal use of the IXFFAtomType class.
  using FFAtom = std::shared_ptr<IXFFAtom>;
  //! \brief shared_ptr for normal use of the IXFFBond class.
  using FFBond = std::shared_ptr<IXFFBond>;
  //! \brief shared_ptr for normal use of the IXFFAngle class.
  using FFAngle = std::shared_ptr<IXFFAngle>;
  //! \brief shared_ptr for normal use of the IXFFDihedral class.
  using FFDihedral = std::shared_ptr<IXFFDihedral>;
  
  /*! \brief Enum for the different families of forcefields.
   *  \details Though other values are provided, currently only the GROMOS
   *  family of forcefields are technically supported. */
  enum class FFFamily {
    Empty,
    GROMOS,
    CHARMM,
    AMBER,
    Other
  };
  
  /*! \brief Enum for the different function types a bond term can have.
   *  \details Though other values are provided, currently only the harmonic and
   *  quartic function types are technically supported. */
  enum class BondType {
    Empty,
    Harmonic,
    Quartic,
    Morse,
    Cubic,
    FENE
  };
  
  /*! \brief Enum for the different function types a angle term can have.
   *  \details Though other values are provided, currently only the harmonic and
   *  cosine-harmonic function types are technically supported. */
  enum class AngleType {
    Empty,
    Harmonic,
    CosineHarmonic,
    UreyBradley,
    Quartic
  };
  
  /*! \brief Enum for the different function types a dihedral term can have.
   *  \details Though other values are provided, currently only the proper and
   *  improper function types are technically supported. */
  enum class DihedralType {
    Empty,
    Proper,
    Improper,
    RyckaertBellemans,
    PeriodicImproper,
    Fourier,
    Restricted
  };
  
  class IXFFDihedral {
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
    /*! \brief Normal construtor
     *  \details Creates a new dihedral term of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function. */
    IXFFDihedral(DihedralType type, int_ id,
                     std::initializer_list<float_> parameters);
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFDihedral>& construct,
                                   const uint32_t version);
    
  public:
    IXFFDihedral() = delete; // no default constructor
    
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
    inline DihedralType GetType() const { return _type; }
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    inline int_ GetID() const { return _id; }
    
  private:
    //! The type of the potential energy function.
    DihedralType _type;
    //! An id for this dihedral type.
    int_ _id;
    //! Terms for the potential energy function.
    DataStore _dat;
    //! Mask for the allowed parameters.
    AllowedMask _mask;
  };
  
  class IXFFAngle {
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
    /*! \brief Normal constructor
     *  \details Creates a new angle of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function. */
    IXFFAngle(AngleType type, int_ id,
                  std::initializer_list<float_> parameters);
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFAngle>& construct,
                                   const uint32_t version);
    
  public:
    IXFFAngle() = delete;   // no default constructor
    
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
    inline AngleType GetType() const { return _type; }
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    inline int_ GetID() const { return _id; }
    
    /*! \brief Get the linked type for this type.
     *  \details A linked type is an equivalent type with a different functional
     *  type.
     *  \return the linked type. */
    inline FFAngle GetLinkedType() const { return _link; }
    
  private:
    //! The type of the potential energy function.
    AngleType _type;
    //! An id for this angle type.
    int_ _id;
    //! Terms for the potential energy function.
    DataStore _dat;
    //! Mask for the allowed parameters.
    AllowedMask _mask;
    //! Linked angle type
    FFAngle _link;
  };
  
  class IXFFBond {
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
    /*! \brief Normal constructor
     *  \details Creates a new angle of the specified type.
     *  \param type the type of the potential energy function to use.
     *  \param id the ID of the type.
     *  \param parameters the list of parameters for the function. */
    IXFFBond(BondType type, int_ id,
                std::initializer_list<float_> parameters);
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFBond>& construct,
                                   const uint32_t version);
    
  public:
    IXFFBond() = delete;   // no default constructor
    
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
    inline BondType GetType() const { return _type; }
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this dihedral type. */
    inline int_ GetID() const { return _id; }
    
    /*! \brief Get the linked type for this type.
     *  \details A linked type is an equivalent type with a different functional
     *  type.
     *  \return the linked type. */
    inline FFBond GetLinkedType() const { return _link; }
    
  private:
    //! The type of the potential energy function.
    BondType _type;
    //! An id for this angle type.
    int_ _id;
    //! Terms for the potential energy function.
    DataStore _dat;
    //! Mask for the allowed parameters.
    AllowedMask _mask;
    //! Linked type
    FFBond _link;
  };
  
  class IXFFAtom {
    //! \brief Friendship allows IXFFAtomType internals to be tested.
    friend struct indigox::test::TestFFAtom;
    //! \brief Friendship allows serialisation
    friend class cereal::access;
    //! \brief Friendship allows atom types to be added to a forcefield.
    friend class IXForcefield;
    
  private:
    /*! \brief Normal constructor
     *  \details Creates a new atom term.
     *  \param id the ID of the type.
     *  \param name the name of the atom type. */
    IXFFAtom(int_ id, string_ name);
    
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXFFAtom>& construct,
                                   const uint32_t version);
    
  public:
    IXFFAtom() = delete;  // no default constructor
    
    /*! \brief Get the ID for this type.
     *  \return the ID of this atom type. */
    inline int_ GetID() const { return _id; }
    
    /*! \brief Get the name for this type
     *  \return the name of this atom type. */
    inline string_ GetName() const { return _name; }
    
  private:
    //! ID of this atom type
    int_ _id;
    //! Name of this atom type
    string_ _name;
  };
  
  class IXForcefield {
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
    IXForcefield(FFFamily family, string_ name);
    
  private:
    template <typename Archive>
    void save(Archive& archive, const uint32_t version) const;
    
    template <typename Archive>
    static void load_and_construct(Archive& archive,
                                   cereal::construct<IXForcefield>& construct,
                                   const uint32_t version);
    
    /*! \brief Adds a new bond type to the forcefield.
     *  \param type the type of the bond.
     *  \param id the id of the bond.
     *  \param parameters the parameters of the bond.
     *  \return the newly created bond type. */
    FFBond NewBondType(BondType type, int_ id,
                       std::initializer_list<float_> parameters);
    
    /*! \brief Adds a new angle type to the forcefield.
     *  \param type the type of the angle.
     *  \param id the id of the angle.
     *  \param parameters the parameters of the angle.
     *  \return the newly created angle type. */
    FFAngle NewAngleType(AngleType type, int_ id,
                             std::initializer_list<float_> parameters);
    
    /*! \brief Adds a new dihedral type to the forcefield.
     *  \param type the type of the dihedral.
     *  \param id the id of the dihedral.
     *  \param parameters the parameters of the dihedral.
     *  \return the newly created dihedral type. */
    FFDihedral NewDihedralType(DihedralType type, int_ id,
                                   std::initializer_list<float_> parameters);
    
  public:
    /*! \brief Reserve space for atom types.
     *  \param size the amount of space to reserve. */
    void ReserveAtomTypes(size_ size) { _atms.reserve(size); }
    
    /*! \brief Reserve space for bond types.
     *  \param type the type of bond to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveBondTypes(BondType type, size_ size) {
      if (_bnds.find(type) != _bnds.end()) _bnds.at(type).reserve(size);
    }
    
    /*! \brief Reserve space for angle types.
     *  \param type the type of angle to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveAngleTypes(AngleType type, size_ size) {
      if (_angs.find(type) != _angs.end()) _angs.at(type).reserve(size);
    }
    
    /*! \brief Reserve space for dihedral types.
     *  \param type the type of dihedral to reserve space for.
     *  \param size the amount of space to reserve. */
    void ReserveDihedralTypes(DihedralType type, size_ size) {
      if (_dhds.find(type) != _dhds.end()) _dhds.at(type).reserve(size);
    }
    
    /*! \brief Get the family of the forcefield
     *  \return the forcefields family. */
    inline FFFamily GetFamily() const { return _family; }
    
    /*! \brief Get the name of the forcefield
     *  \return the name of the forcefield. */
    inline string_ GetName() const { return _name; }
    
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
    
    /*! \brief Add a new harmonic bond type.
     *  \param id the id of the type.
     *  \param k the force constant of the type.
     *  \param length the ideal bond length of the type.
     *  \return the newly created bond type. */
    FFBond NewHarmonicBondType(int_ id, float_ k, float_ length) {
      return NewBondType(BondType::Harmonic, id, {k, length});
    }
    
    /*! \brief Add a new quartic bond type.
     *  \param id the id of the type.
     *  \param k the force constant of the type.
     *  \param length the ideal bond length of the type.
     *  \return the newly created bond type. */
    FFBond NewQuarticBondType(int_ id, float_ k, float_ length) {
      return NewBondType(BondType::Quartic, id, {k, length});
    }
    
    /*! \brief Get the type of bond eith the given id.
     *  \param type the type of bond to get.
     *  \param the id of the bond type.
     *  \return the requested bond type, or an empty shared_ptr. */
    FFBond GetBondType(BondType type, int_ id) const;
    
    /*! \brief Get the harmonic bond type with the given id.
     *  \param id the id of the bond type.
     *  \return the harmonic bond type, or an empty shared_ptr. */
    FFBond GetHarmonicBondType(int_ id) const {
      return GetBondType(BondType::Harmonic, id);
    }
    
    /*! \brief Get the quartic bond type with the given id.
     *  \param id the id of the bond type.
     *  \return the quartic bond type, or an empty shared_ptr. */
    FFBond GetQuarticBondType(int_ id) const {
      return GetBondType(BondType::Quartic, id);
    }
    
    /*! \brief Add a new harmonic angle type.
     *  \param id the id of the type.
     *  \param k the force constant of the type.
     *  \param theta the ideal bond angle of the type.
     *  \return the newly created angle type, or an empty shared_ptr. */
    FFAngle NewHarmonicAngleType(int_ id, float_ k, float_ theta) {
      return NewAngleType(AngleType::Harmonic, id, {k, theta});
    }
    
    /*! \brief Link together two equivalent bond types.
     *  \param a,b the two types to link. */
    void LinkBondTypes(FFBond a, FFBond b) {
      a->_link = b; b->_link = a;
    }
    
    /*! \brief Link together two equivalent angle types.
     *  \param a,b the two types to link. */
    void LinkAngleTypes(FFAngle a, FFAngle b) {
      a->_link = b; b->_link = a;
    }
    
    /*! \brief Add a new cosine-harmonic angle type.
     *  \param id the id of the type.
     *  \param k the force constant of the type.
     *  \param theta the ideal bond angle of the type.
     *  \return the newly created angle type, or an empty shared_ptr. */
    FFAngle NewCosineHarmonicAngleType(int_ id, float_ k, float_ theta) {
      return NewAngleType(AngleType::CosineHarmonic, id, {k, theta});
    }
    
    /*! \brief Get the type of angle with the given id.
     *  \param type the type of angle to get.
     *  \param the id of the angle type.
     *  \return the requested angle type, or an empty shared_ptr. */
    FFAngle GetAngleType(AngleType type, int_ id) const;
    
    /*! \brief Get the harmonic angle type with the given id.
     *  \param id the id of the angle type.
     *  \return the harmonic angle type, or an empty shared_ptr. */
    FFAngle GetHarmonicAngleType(int_ id) const {
      return GetAngleType(AngleType::Harmonic, id);
    }
    
    /*! \brief Get the cosine-harmonic angle type with the given id.
     *  \param if the id of the angle type.
     *  \return the cosine-harmonic angle type or an empty shared_ptr. */
    FFAngle GetCosineHarmonicAngleType(int_ id) const {
      return GetAngleType(AngleType::CosineHarmonic, id);
    }
    
    /*! \brief Add a new harmonic angle type.
     *  \param id the id of the type.
     *  \param k the force constant of the type.
     *  \param phase the phase shift of the type.
     *  \param m the multiplicity of the type.
     *  \return the newly created angle type, or an empty shared_ptr. */
    FFDihedral NewProperDihedralType(int_ id, float_ k, float_ phase, uint_ m) {
      return NewDihedralType(DihedralType::Proper,
                             id, {phase, k, static_cast<float_>(m)});
    }
    
    /*! \brief Add a new improper dihedral type.
     *  \param id the id of the type.
     *  \param k the force constant of the type.
     *  \param theta the ideal angle of the type.
     *  \return the newly created dihedral type, or an empty shared_ptr. */
    FFDihedral NewImproperDihedralType(int_ id, float_ k, float_ theta) {
      return NewDihedralType(DihedralType::Improper, id, {theta, k});
    }
    
    /*! \brief Get the type of dihedral with the given id.
     *  \param type the type of dihedral to get.
     *  \param the id of the dihedral type.
     *  \return the requested dihedral type, or an empty shared_ptr. */
    FFDihedral GetDihedralType(DihedralType type, int_ id) const;
    
    /*! \brief Get the improper dihedral type with the given id.
     *  \param id the id of the dihedral type.
     *  \return the improper dihedral type, or an empty shared_ptr. */
    FFDihedral GetImproperDihedralType(int_ id) const {
      return GetDihedralType(DihedralType::Improper, id);
    }
    
    /*! \brief Get the proper dihedral type with the given id.
     *  \param id the id of the dihedral type.
     *  \return the proper dihedral type, or an empty shared_ptr. */
    FFDihedral GetProperDihedralType(int_ id) const {
      return GetDihedralType(DihedralType::Proper, id);
    }
    
    /*! \brief Get the number of atom types in the forcefield.
     *  \return the number of atom types in the forcefield. */
    size_ NumAtomTypes() const { return _atms.size(); }
    
    /*! \brief Get the number of bond types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of bond types in the forcefield. */
    size_ NumBondTypes() const {
      return std::accumulate(_bnds.begin(), _bnds.end(), 0,
                             [](auto a, auto& b){return a + b.second.size();});
    }
    
    /*! \brief Get the number of angle types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of angle types in the forcefield. */
    size_ NumAngleTypes() const {
      return std::accumulate(_angs.begin(), _angs.end(), 0,
                             [](auto a, auto& b){return a + b.second.size();});
    }
    
    /*! \brief Get the number of dihedral types in the forcefield.
     *  \details For forcefields with multiple allowed types, the count of each
     *  type are summed.
     *  \return the number of bond types in the forcefield. */
    size_ NumDihedralTypes() const {
      return std::accumulate(_dhds.begin(), _dhds.end(), 0,
                             [](auto a, auto& b){return a + b.second.size();});
    }
    
  private:
    //! The family of the forcefield
    FFFamily _family;
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
  
  Forcefield GenerateGROMOS54A7();
}

#endif    /* INDIGOX_CLASSES_FORCEFIELD_HPP */
