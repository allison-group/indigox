/*! \file periodictable.hpp
 */

#ifndef INDIGOX_CLASSES_PERIODIC_TABLE_HPP
#define INDIGOX_CLASSES_PERIODIC_TABLE_HPP

#include <cstdint>
#include <iostream>
#include <map>
#include <string>

#include "../utils/common.hpp"

namespace indigox {
  
  class IXPeriodicTable;
  class IXElement;
  
  
  //! \brief shared_ptr for normal use of the IXPeriodicTable class.
  typedef std::shared_ptr<IXPeriodicTable> PeriodicTable;
  typedef std::shared_ptr<IXElement> Element;
  /*! \brief weak_ptr for non-ownership reference to the IXPeriodicTable class.
   *  \details Intended for internal use only. */
  typedef std::weak_ptr<IXPeriodicTable> _PeriodicTable;
  typedef std::weak_ptr<IXElement> _Element;
  
  /*! \class IXPeriodicTable periodictable.hpp indigox/classes/periodictable.hpp
   *  \brief Singleton class for storing and access elemental information.
   *  \details The IXPeriodicTable class provides the only means to access
   *  the IXElement class. Access to the instance should only be obtained
   *  using the GetInstance() method. Like most other classes in the indigoX
   *  library, usage is primarily intended through the use of smart pointers. */
  class IXPeriodicTable {
  public:
    
    /*! \brief Obtain the singleton instance of the PeriodicTable.
     *  \details If an instance does not exist, creates one and generates
     *  elemental information.
     *  \returns the PeriodicTable instance. */
    static PeriodicTable GetInstance();
    
  public:
    /*! \brief Get the element with the given atomic number.
     *  \param z the atomic number of the element to get.
     *  \returns the requested element.
     *  \throw std::invalid_argument If the requested atomic number does not exist
     *  within the PeriodicTable. */
    Element GetElement(const uint8_t z) const;
    
    /*! \brief Get the element with the given name or symbol.
     *  \details Matches are made ignoring case.
     *  \param name the name or symbol of the element to get.
     *  \returns the requested element.
     *  \throw std::invalid_argument If the requested name or symbol does not
     *  exist within the PeriodicTable. */
    Element GetElement(const std::string name) const;
    
    /*! \brief Get the element with the given atomic number.
     *  \param z the atomic number of the element to get.
     *  \return the requested element.
     *  \see IXPeriodicTable::GetElement */
    Element operator[](const uint8_t z) const { return GetElement(z); }
    
    /*! \brief Get the element with the given name or symbol.
     *  \param name the name or symbol of the element to get.
     *  \return the requested element.
     *  \see IXPeriodicTable::GetElement */
    Element operator[](const std::string name) const { return GetElement(name); }
    
    /*! \brief Get the element for use when an element is not defined.
     *  \details As there is not much point in an undefined element, this
     *  method is intended for internal use.
     *  \return the undefined Element. */
    Element GetUndefinedElement() const { return _null; }
    
    /*! \brief Number of elements in the PeriodicTable.
     *  \return the number of elements in the PeriodicTable. */
    size_t NumElements() const { return _z_to.size(); }
    
    /*! \brief Get a textual representation of the PeriodicTable.
     *  \details The textual representation lays out the elements of the
     *  PeriodicTable as one would expect in a normal periodic table. Each
     *  element has its atomic number and symbol printed.
     *  \return a string containing a tabular layout of the PeriodicTable
     *  elements. */
    std::string ToString() const;
  
  private:
    //! Default constructor.
    IXPeriodicTable() = default;
    
    /*! \brief Generate Elements.
     *  \details Currently, the first 109 elements are generated. These are
     *  hard coded into the implementation file. These are populated into both
     *  the _z_to and _name_to maps. */
    void GeneratePeriodicTable();
  
  private:
    //! Has been initalised.
    static bool _init;
    //! The initailied instance.
    static PeriodicTable _instance;
    //! Undefined element.
    Element _null;
    //! Map atomic numbers to elements.
    std::map<uint8_t, Element> _z_to;
    //! Map symbol and name to elements.
    std::map<std::string, Element> _name_to;
  };
  
  /** @class IXElement periodictable.hpp classes/periodictable.hpp
   *  @brief Read only class for storing elemental information.
   *  @details Contains a large amount of relevant information pertaining
   *  to elements. No public constructors so can only be created by the
   *  PeriodicTable class.
   *  @since 0.1
   */
  class IXElement {
  public:
    /// @returns the relative atomic mass of the element in daltons.
    double GetAtomicMass() const { return mass_; }
    
    /// @returns the atomic number of the element.
    uint8_t GetAtomicNumber() const { return Z_; }
    
    /// @returns the atomic radius of the element in angstroms.
    double GetAtomicRadius() const { return radius_; }
    
    /// @returns the covalent radius of the element in angstroms.
    double GetCovalentRadius() const { return cov_; }
    
    /// @returns the van der Waals radius of the element in angstroms.
    double GetVanDerWaalsRadius() const { return vdw_; }
    
    /// @returns the name of the element.
    std::string GetName() const { return name_; }
    
    /// @returns the symbol of the element.
    std::string GetSymbol() const { return symbol_; }
    
    /// @returns the IUPAC group number the element is in.
    uint8_t GetGroup() const { return grp_; }
    
    /// @returns the period the element is in.
    uint8_t GetPeriod() const { return period_; }
    
    /// @returns the number of valence electrons the element contains.
    uint8_t GetValenceElectronCount() const { return val_; }
    
    /// @returns the number of electrons required for a full outer shell.
    uint8_t GetOctet() const { return oct_; }
    
    /// @returns the number of electrons allowed when allowing hypervalency.
    uint8_t GetHypervalentOctet() const { return hyper_; }
    
    /// @returns the electronegativity of the element on the Pauling scale.
    double GetElectronegativity() const { return chi_; }

    /// @returns a simple string representation of the element.
    std::string ToString() const;
    
  private:
    const std::string name_, symbol_;
    const uint8_t grp_, period_, Z_, val_, oct_, hyper_;
    const double mass_, radius_, cov_, vdw_, chi_;
    
  private:  // Only create Elements within PeriodicTable
    friend class IXPeriodicTable;
    IXElement();
    IXElement(uint8_t, std::string, std::string, double, uint8_t, uint8_t,
              uint8_t, uint8_t, uint8_t, double, double, double, double);
  };
} // namespace indigox

// Operators have to be explicitly inlined or python bindings linkage fails
inline std::ostream& operator<<(std::ostream& os, indigox::Element e) {
  return e ? (os << e->ToString()) : os;
}

/*! \brief Print a PeriodicTable to an output stream.
 *  \details Prints number of elements in the PeriodicTable.
 *  \param os output stream to print to.
 *  \param pt PeriodicTable to print.
 *  \return the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, indigox::PeriodicTable pt) {
  return pt ? (os << "PeriodicTable(" << pt->NumElements() << " elements)") : os;
}

#endif /* INDIGOX_CLASSES_PERIODIC_TABLE_HPP */
