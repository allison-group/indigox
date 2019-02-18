#include "../classes/athenaeum.hpp"
#include "../classes/forcefield.hpp"
#include "../utils/enum_class_bitwise.hpp"
#include "../utils/fwd_declares.hpp"

#include <bitset>
#include <cstdint>
#include <list>
#include <type_traits>

#ifndef INDIGOX_ALGORITHM_CHERRYPICKER_HPP
#define INDIGOX_ALGORITHM_CHERRYPICKER_HPP

namespace indigox::algorithm {

  /*! \brief CherryPicker parameterisation algorithm class.

   The CherryPicker algorithm is a parameterisation algorithm for molecular
   dynamics simulations. It parameterises novel molecules by identifying
   fragments of previously parameterised molecules which match portions of the
   novel molecule, then applying the parameters found to the novel molecule. See
   the \ref usingcherrypicker "tutorial" for a brief tutorial on how to
   parameterise something using the CherryPicker algorithm.

   There are three data classes and one graph class primarily associated with
   CherryPicker. The Athenaeum class acts as a library of previously
   parameterised molecules which can be used to parameterise a target molecule.
   Within an Athenaeum, each Molecule is mapped to a vector of Fragment
   instances. The Fragment class represents a fragment of a source molecule. It
   is the fragments which are checked for matching portions of the target
   molecule. This checking is performed using the \link
   graph::CondensedMolecularGraph CondensedMolecularGraph (CMG)\endlink of the
   target molecule and the fragment. Accounting of the matches found is
   performed by the ParamMolecule class and its components.

   To parameterise a novel molecule using CherryPicker, first an Athenaeum must
   be created (see the Athenaeum documentation). For a given Athenaeum, the
   process of parameterisation is as follows. Fragments within the Athenaeum are
   iterated through in an arbitary order. The CMG of each fragment is checked
   against the CMG of the target molecule for subgraph isomorphism. All subgraph
   isomorphisms between the two graphs are enumerated, giving a series of
   one-to-one mappings between the vertices of each graph. These mappings are
   expanded to give one-to-one mappings between the atoms of the fragment's
   source molecule and the target molecule. The use of all possible mappings
   between the condensed atoms of each vertex is available through the Settings.
   From these mappings, the atoms, bonds, angles and dihedrals that match one
   another from the fragment source molecule and the target molecule are known.
   Parameters assigned to these components in the fragment source molecule are
   accounted for by adding them to the corresponding component within the
   ParamMolecule created from the target molecule. Once all fragments within an
   Athenaeum have been checked, any component within the ParamMolecule which has been mapped by at least one component from a fragment applies a parameter back to the corresponding component within the target molecule. In the case of discrete parameters, such as bond types, the parameter applied is the mode of all mapped parameters. In the case of continuous parameters, such as atomic partial charges, the parameter applied is the mean of all mapped parameters. If an Athenaeum is marked as being self-consistent, it is at the point of applying these parameters that the self-consistency is checked. For discrete parameters, self-consistent is defined as having only one mapped type on the ParamMolecule component. For continuous parameters, it is defined as having a difference of no greater than \f$10^{-10}\f$ between the largest and smallest mapped values.

   CherryPicker is designed to support multiple Athenaeums. This is handled in a
   first in-first out (FIFO). That is, Athenaeums are iterated through in the
   order in which they were assigned to the CherryPicker instance. When multiple Athenaeums are provided, parameterisation occurs as above for each Athenaeum. Additionally, once all fragments within an Athenaeum have been checked and parameters applied, any components which had a parameter applied is removed from the set of parameterisable components so that no further Athenaeums will be able to parameterise it. The reasoning for this FIFO methadology is so that portions of a molecule can be parameterised by a smaller set of well tested fragments. For example, a protein with a few unnatural or unique residues can initially be parameterised with fragments matching the natural amino acids, thus covering the majority of the molecule quickly, before a more shotgun approach to parameterising the remaining portions is utilised.
   
   An exception to the rule of not further parameterising a component which has been parameterised by a previous Athenaeum is when any excess charge is to redistributed across the molecule. In this case, all atoms which were not parameterised by a self-consistent Athenaeum are available to have extra charge added.

   The primary working component of the algorithm is <a
   href="https://en.wikipedia.org/wiki/Subgraph_isomorphism_problem">subgraph
   isomorphism</a>. Whenever a fragment is tested to check if it matches the
   test molecule, the CMG of the fragment is tested for subgraph isomorphism
   with the CMG of the target molecule.

    There are two kinds of settings, boolean and integer values. Attempting to
   modify a setting through the incorrect access type (line 4 of the above code
   snippet) will result in an error being thrown. See the \link Settings
   Settings\endlink documentation for details.

    The ParameteriseMolecule() method applies the CherryPicker algorithm to the
   provided \link Molecule molecule\endlink. Though it returns a ParamMolecule,
   it is safe to ignore this as the parameters discovered are applied directly
   to the provided \link Molecule molecule\endlink. If you desire to modify the
   parameters discovered, the returned ParamMolecule provides all of the
   information relating to what parameters were found for each component of the
   \link Molecule molecule\endlink. See the Molecule class for serialisation
   options.
   */
  class CherryPicker {
  public:
    /*! \brief User controllable settings for the CherryPicker algorithm.

    Settings are maintained on a per-instance basis. That is multiple different
    CherryPicker instances can have different settings. There are two types of
    settings, boolean and integer. Boolean settings are either on or off, and
     are manipulated through the SetBool() and UnsetBool() methods. These
    settings appear before the \link Settings::BoolCount BoolCount\endlink
    marker. Integer settings require an integer value and are manipulated
    through the SetInt() method. These settings appear between the \link
    Settings::BoolCount BoolCount\endlink and \link Settings::IntCount
    IntCount\endlink markers. Default values are detailed in the
    DefaultSettings() method.
     */
    enum class Settings : uint8_t {
      // Vertex
      /*! When performing subgraph isomorphism testing, two vertices will only
         match if the element of their associated atoms are the same. */
      VertexElement,
      /*! When performing subgraph isomorphism testing, two vertices will only
         match if the formal charge of their associated atoms are the same. */
      VertexFormalCharge,
      /*! When performing subgraph isomorphism testing, two vertices will only
         match if the counts per type of condensed vertices are the same. */
      VertexCondensed,
      /*! When performing subgraph isomorphism testing, two vertices will only
       match if they are both contained within a cycle. To be regarded as in a
       cycle, the size of the smallest cycle a vertex is in must be less than or
       equal to the ??? setting.

       \todo Make the setting in CMG. */
      VertexCyclic,
      /*! When performing subgraph isomorphism testing, two  vertices will only
         match if the smallest cycle they are both contained in have the same
         size. */
      VertexCyclicSize,
      /*! When performing subgraph isomorphism testing, two vertices will only
         match if the assocaited atoms have the same stereochemistry. */
      VertexStereochemistry,
      /*! When performing subgraph isomorphism testing, two vertices will only
         match if the assocaited atoms have the same aromaaticity. */
      VertexAromaticity,
      /*! When performing subgraph isomorphism testing, two vertices will only
         match if the assocaited atoms have the same number of bonds. */
      VertexDegree,
      // Edge
      /*! When performing subgraph isomorphism testing, two \link graph::CMGEdge
         edges\endlink will only match if the \link Bond::Order bond
         order\endlink of their associated \link Bond bonds\endlink are the
         same.*/
      EdgeBondOrder,
      /*! When performing subgraph isomorphism testing, two \link graph::CMGEdge
         edges\endlink will only match if their assocaited \link Bond
         bonds\endlink have the same \link Bond::Stereo stereochemistry\endlink.
       */
      EdgeStereochemistry,
      /*! When performing subgraph isomorphism testing, two \link graph::CMGEdge
       edges\endlink will only match if they are both contained within a cycle.
       To be regarded as in a cycle, the size of the smallest cycle a \link
       graph::CMGEdge edge\endlink is in must be less than or equal to the ???
       setting. \todo Make the setting in CMG. */
      EdgeCyclic,
      /*! When performing subgraph isomorphism testing, two \link graph::CMGEdge
         edges\endlink will only match if the smallest cycle they are both
         contained in have the same size. */
      EdgeCyclicSize,
      /*! When performing subgraph isomorphism testing, two \link graph::CMGEdge
         edges\endlink will only match if the \link Atom atoms\endlink of the
         associated \link Bond bond\endlink have the same number of \link Bond
         bonds\endlink. That is, the \link Atom atom\endlink of each \link Bond
         bond\endlink with the lowest number of \link Bond bonds\endlink must
         have the same number of \link Bond bonds\endlink, and the same for the
         \link Atom atom\endlink with the highest number of \link Bond
         bonds\endlink. */
      EdgeDegree,
      // Normal bool
      /*! When applying parameters to a \link ParamBond bond\endlink, allow for
         parameteristion across the overlap region. If one \link Atom
         atom\endlink of a \link Bond bond\endlink is in the \link Fragment
         fragment\endlink region and the other is in the overlap, setting this
         allows the matching \link ParamBond bond\endlink to be parameterised by
         the \link Bond bond\endlink between the fragment and overlap regions.
       */
      AllowDanglingBonds,
      /*! When applying parameters to an \link ParamAngle angle\endlink, allow
         for parameteristion across the overlap region. If two \link Atom
         atoms\endlink of a \link Angle angle\endlink are in the \link Fragment
         fragment\endlink region and the other is in the overlap, setting this
         allows the matching \link ParamAngle angle\endlink to be parameterised
         by the \link Angle angle\endlink between the fragment and overlap
         regions. */
      AllowDanglingAngles,
      /*! When applying parameters to a \link ParamDihedral dihedral\endlink,
         allow for parameterisation across the overlap region. If two or three
         \link Atom atoms\endlink of a \link Dihedral dihedral\endlink are in
         the \link Fragment fragment\endlink region and the others are in the
         overlap region, setting this allows the matching \link ParamDihedral
         dihedral\endlink to be parameterised by the \link Dihedral
         dihedral\endlink between the fragment and overlap regions. This only
         applies to \link Dihedral dihedrals\endlink where the \link Atom
         atoms\endlink that are in the fragment and not the overlap are adjacent
         to one another. */
      AllowDanglingDihedrals,
      /*! When applying parameters iterate through all permutations of
       contracted vertex mappings. On any given \link graph::CMGVertex
       vertex\endlink match, \link graph::MGVertex contracted vertices\endlink
       of the same \link graph::CMGVertex::ContractedSymmetry type\endlink can
       be mapped to one another in all possible permutations. This setting
       causes all such permutations to be used for applying parameters. In
       general, this is not required as \link graph::MGVertex contracted
       vertices\endlink are expected to have the same parameters.
       \note Enabling this option will cause computational requirements to
       increase dramatically. */
      ParameteriseFromAllPermutations,
      /*! Marks the end of the boolean settings. As there is no external use for
         this value, it is not exposed to Python. */
      BoolCount,
      /*! The minimum size of \link Fragment fragments\endlink to check for
         matching purposes. Any \link Fragment fragment\endlink within an \link
         Athenaeum athenaeum\endlink which has a size smaller than this will be
         skipped. A negative value means all \link Fragment fragments\endlink
         will be tested. */
      MinimumFragmentSize,
      /*! The maximum size of fragments to check for matching purposes. Any
         fragment in an athenaeum which has a size larger than this will be
         skipped. A negative value means all fragments will be tested. */
      MaximumFragmentSize,
      /*! Marks the end of the integer settings. As there is no external use for
         this value, it is not exposed to Python. */
      IntCount
    };

    /*! \brief Get the current state of a boolean setting.

     \param param The setting to get the state of.
     \returns the current state of the \p param setting.
     \throws std::runtime_error if \p param is not a valid boolean setting.
     */
    bool GetBool(Settings param);

    /*! \brief Set the state of a boolean setting to true.

     \param param The setting to set the state of.
     \throws std::runtime_error if \p param is not a valid boolean setting.
     */
    void SetBool(Settings param);

    /*! \brief Set the state of a boolean setting to false.

     \param param The setting to set the state of.
     \throws std::runtime_error if \p param is not a valid boolean setting.
     */
    void UnsetBool(Settings param);

    /*! \brief Get the current value of an integer setting.

     \param param The setting to get the value of.
     \returns the current value of the \p param setting.
     \throws std::runtime_error if \p param is not a valid integer setting.
     */
    int32_t GetInt(Settings param);

    /*! \brief Set the value of an integer setting.

     \param param The setting to set the value of.
     \param value The value to set \p param to.
     \throws std::runtime_error if \p param is not a valid integer setting.
     */
    void SetInt(Settings param, int32_t value);

    /*! \brief Set the default values for all of the settings.

     The following boolean values are default set to true:
     \link Settings::VertexElement VertexElement\endlink, \link
     Settings::VertexFormalCharge VertexFormalCharge\endlink, \link
     Settings::VertexCondensed VertexCondensed\endlink, \link
     Settings::VertexDegree VertexDegree\endlink, \link Settings::EdgeBondOrder
     EdgeBondOrder\endlink, \link Settings::EdgeDegree EdgeDegree\endlink, \link
     Settings::AllowDanglingBonds AllowDanglingBonds\endlink, \link
     Settings::AllowDanglingAngles AllowDanglingAngles\endlink, and \link
     Settings::AllowDanglingDihedrals AllowDanglingDihedrals\endlink. All other
     boolean values default to false. The default \link
     Settings::MinimumFragmentSize MinimumFragmentSize\endlink is \f$4\f$ and
     the default \link Settings::MaximumFragmentSize MaximumFragmentSize\endlink
     is \f$-1\f$.
     */
    void DefaultSettings();

    /*! \brief Constructor

     \param ff defines what forcefield will be used for parameterisation.
     */
    CherryPicker(Forcefield &ff);

    /*! \brief Add an Athenaeum for parameterisation purposes.

     Adds an Athenaeum in a FIFO manner. No check is performed to see if \p
     library has previously been added, so multiple instances of the same
     Athenaeum is possible. Only check performed is that the Athenaeum's
     forcefield matches the CherryPicker forcefield.

     \param library the Athenaeum to add.
     \returns if \p library was successfully added.
     */
    bool AddAthenaeum(Athenaeum &library);

    /*! \brief Remove an Athenaeum from the list.

     \param library the Athenaeum to remove.
     \returns if \p library was successfully removed.
     */
    bool RemoveAthenaeum(Athenaeum &library);

    /*! \brief The number of Athenaeums in the list.

     \returns the number of Athenaeums.
     */
    int32_t NumAthenaeums() { return (int32_t)_libs.size(); }

    /*! \brief Apply the CherryPicker algorithm to a molecule.

     Fill in the details here.

     \param mol the molecule to parameterise.
     \returns a ParamMolecule for all the matched parameters.
     \throws std::runtime_error if the list of Athenaeums is empty or \p mol is
     not connected.
     */
    ParamMolecule ParameteriseMolecule(Molecule &mol);

    /*! \brief Get the assigned forcefield.

     \returns the assigned Forcefield.
     */
    Forcefield GetForcefield() { return _ff; }

    //! \cond ignored
    CherryPicker() = delete;
    ~CherryPicker() = default;
    //! \endcond

  private:
    Forcefield _ff;
    std::list<Athenaeum> _libs;
    std::bitset<(uint8_t)Settings::BoolCount> bool_parameters;
    std::array<int32_t,
               (uint8_t)Settings::IntCount - (uint8_t)Settings::BoolCount - 1>
        int_parameters;
  };

} // namespace indigox::algorithm

#endif /* INDIGOX_ALGORITHM_CHERRYPICKER_HPP */
