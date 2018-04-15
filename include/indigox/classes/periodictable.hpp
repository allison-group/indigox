/*! \file periodictable.hpp
 */

#ifndef INDIGOX_CLASSES_PERIODIC_TABLE_HPP
#define INDIGOX_CLASSES_PERIODIC_TABLE_HPP

#include <iostream>
#include <map>

#include "../utils/common.hpp"
#include "../utils/numerics.hpp"

namespace indigox {
  
  class IXPeriodicTable;
  class IXElement;
  
  
  //! \brief shared_ptr for normal use of the IXPeriodicTable class.
  typedef std::shared_ptr<IXPeriodicTable> PeriodicTable;
  //! \brief shared_ptr for normal use of the IXElement class.
  typedef std::shared_ptr<IXElement> Element;
  /*! \brief weak_ptr for non-ownership reference to the IXPeriodicTable class.
   *  \details Intended for internal use only. */
  typedef std::weak_ptr<IXPeriodicTable> _PeriodicTable;
  /*! \brief weak_ptr for non-ownership reference to the IXElement class.
   *  \details Intended for internal use only. */
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
    Element GetElement(const uchar_ z) const;
    
    /*! \brief Get the element with the given name or symbol.
     *  \details Matches are made ignoring case.
     *  \param name the name or symbol of the element to get.
     *  \returns the requested element.
     *  \throw std::invalid_argument If the requested name or symbol does not
     *  exist within the PeriodicTable. */
    Element GetElement(const string_ name) const;
    
    /*! \brief Get the element with the given atomic number.
     *  \param z the atomic number of the element to get.
     *  \return the requested element.
     *  \see IXPeriodicTable::GetElement(const uint8_t) const */
    Element operator[](const uchar_ z) const { return GetElement(z); }
    
    /*! \brief Get the element with the given name or symbol.
     *  \param name the name or symbol of the element to get.
     *  \return the requested element.
     *  \see IXPeriodicTable::GetElement(const std::string) const */
    Element operator[](const string_ name) const { return GetElement(name); }
    
    /*! \brief Get the element for use when an element is not defined.
     *  \details As there is not much point in an undefined element, this
     *  method is intended for internal use.
     *  \return the undefined Element. */
    Element GetUndefinedElement() const { return _null; }
    
    /*! \brief Number of elements in the PeriodicTable.
     *  \return the number of elements in the PeriodicTable. */
    size_ NumElements() const { return _z_to.size(); }
    
    /*! \brief Get a textual representation of the PeriodicTable.
     *  \details The textual representation lays out the elements of the
     *  PeriodicTable as one would expect in a normal periodic table. Each
     *  element has its atomic number and symbol printed.
     *  \return a string containing a tabular layout of the PeriodicTable
     *  elements. */
    string_ ToString() const;
  
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
    std::map<uchar_, Element> _z_to;
    //! Map symbol and name to elements.
    std::map<string_, Element> _name_to;
  };
  
  /*! \class IXElement periodictable.hpp indigox/classes/periodictable.hpp
   *  \brief Read only class for storing elemental information.
   *  \details Contains a large amount of relevant information pertaining
   *  to elements. Information includes atomic mass, atomic number, atomic
   *  radius, covalent radius, van der Waals radius, name, symbol, IUPAC group
   *  number, period number, number of valence electrons, octet, hypervalency
   *  octet and electronegativity. No public constructors are provided, so
   *  instances can only be generated by the IXPeriodicTable class. Like most
   *  other classes in the indigoX library, usage is primarily intended through
   *  the use of smart pointers.
   */
  class IXElement {
  public:
    
    /*! \brief Get atomic mass.
     *  \return the atomic mass of the element. */
    float_ GetAtomicMass() const { return _mass; }
    
    /*! \brief Get atomic number.
     *  \return the atomic number of the element. */
    uchar_ GetAtomicNumber() const { return _Z; }
    
    /*! \brief Get atomic radius.
     *  \details Radius is in angstroms.
     *  \return the atomic radius of the element. */
    float_ GetAtomicRadius() const { return _rad; }
    
    /*! \brief Get covalent radius.
     *  \details Radius is in angstroms.
     *  \return the covalent radius of the element. */
    float_ GetCovalentRadius() const { return _cov; }
    
    /*! \brief Get van der Waals radius.
     *  \details Radius is in angstroms.
     *  \return the van der Waals radius of the element. */
    float_ GetVanDerWaalsRadius() const { return _vdw; }
    
    /*! \brief Get element name.
     *  \return the name of the element. */
    string_ GetName() const { return _nme; }
    
    /*! \brief Get element symbol.
     *  \return the symbol of the element. */
    string_ GetSymbol() const { return _sym; }
    
    /*! \brief Get element group number.
     *  \return the IUPAC group number of the element. */
    uchar_ GetGroup() const { return _grp; }
    
    /*! \brief Get element period.
     *  \return the period of the periodic table the element is in. */
    uchar_ GetPeriod() const { return _prd; }
    
    /*! \brief Get number of valence electrons
     *  \return the number of valence electrons the element contains. */
    uchar_ GetValenceElectronCount() const { return _val; }
    
    /*! \brief Get full outer shell octet.
     *  \return the number of electrons required for a full outer shell. */
    uchar_ GetOctet() const { return _oct; }
    
    /*! \brief Get full outer shell octet whne allowing for hypervalency.
     *  \return the number of electrons required for a full outer shell in a
     *  hypervalent state. */
    uchar_ GetHypervalentOctet() const { return _hyp; }
    
    /*! \brief Get electronegativity.
     *  \return the electronegativity of the element on the Pauling scale. */
    float_ GetElectronegativity() const { return _chi; }

    /*! \brief Get a textual representation of the element.
     *  \details Representation contains the element's name, symbol and atomic
     *  number.
     *  \return textutal representation of the element. */
    string_ ToString() const;
    
    //! \cond
    operator bool() { return bool(_Z); }
    //! \endcond
    
  private:
    /*! \property _nme
     *  \brief Element name.
     *  \property _sym
     *  \brief Atomic symbol. */
    const string_ _nme, _sym;
    /*! \property _grp
     *  \brief IUPAC group.
     *  \property _prd
     *  \brief Period.
     *  \property _Z
     *  \brief Atomic number.
     *  \property _val
     *  \brief Number of valence electrons.
     *  \property _oct
     *  \brief Octet.
     *  \property _hyp
     *  \brief Hypervalent octet. */
    const uchar_ _grp, _prd, _Z, _val, _oct, _hyp;
    /*! \property _mass
     *  \brief Relative atomic mass.
     *  \property _rad
     *  \brief Atomic radius.
     *  \property _cov
     *  \brief Covalent radius.
     *  \property _vdw
     *  \brief Van der Waals radius.
     *  \property _chi
     *  \brief Electronegativity. */
    const float_ _mass, _rad, _cov, _vdw, _chi;
    
  private:
    //! Allow IXPeriodicTable to create new IXElement instances.
    friend class IXPeriodicTable;
    
    // No default constructor provided
    IXElement() = delete;
    
    /*! \brief Construct new IXElement instance given data.
     *  \param Z atomic number.
     *  \param name element name.
     *  \param sym atomic symbol.
     *  \param mass relative atomic mass.
     *  \param grp IUPAC group number.
     *  \param prd period.
     *  \param val valence electrons.
     *  \param oct octet.
     *  \param hyp hypervalent octet.
     *  \param rad atomic radius in angstroms.
     *  \param cov covalent radius in angstroms.
     *  \param vdw van der Waals radius i angstroms.
     *  \param chi electronegativity. */
    IXElement(uchar_ Z, string_ name, string_ sym, float_ mass,
              uchar_ grp, uchar_ prd, uchar_ val, uchar_ oct, uchar_ hyp,
              float_ rad, float_ cov, float_ vdw, float_ chi);
  };

  // Operators have to be explicitly inlined or python bindings linkage fails
  
  /*! \brief Print an Element to an output stream.
   *  \details Prints only the element name.
   *  \param os output stream to print to.
   *  \param e Element to print.
   *  \return the output stream. */
  inline std::ostream& operator<<(std::ostream& os, Element e) {
    return e ? (os << "Element(" << e->GetName() << ")") : os;
  }
  
  /*! \brief Print a PeriodicTable to an output stream.
   *  \details Prints number of elements in the PeriodicTable.
   *  \param os output stream to print to.
   *  \param pt PeriodicTable to print.
   *  \return the output stream. */
  inline std::ostream& operator<<(std::ostream& os, PeriodicTable pt) {
    return pt ? (os << "PeriodicTable(" << pt->NumElements() << " elements)") : os;
  }
  
  // Comparison operators
  /*! \brief Equality test of Element and integer.
   *  \details Compares if elements atomic number is equal to the integer. If
   *  the Element is an empty pointer or the undefined element, result will
   *  always be false.
   *  \param l, r element and integer to compare.
   *  \return if the element has the given atomic number. */
  inline bool operator==(Element l, uint8_t r) {
    return (l && *l) ? l->GetAtomicNumber() == r : false;
  }
  /*! \brief Equality test of Element and string.
   *  \details If the size of r is less than or equal to two, compares if the
   *  element has the given atomic symbol. This comparison is case-sensitive.
   *  Otherwise the comparison checks if the element has the given name. This
   *  comparison is not case sensitive. If the Element is an empty pointer or
   *  the undefined element, result will always be false.
   *  \param l, r element and string to compare.
   *  \return if the element has the given atomic symbol or name. */
  inline bool operator==(Element l, std::string r) {
    if (!l || !*l) return false;
    else return r.size() <= 2 ? (l->GetSymbol() == r)
      : (l->GetName() == utils::toUpperFirst(&r));
  }
  /*! \brief Equality test of two Elements.
   *  \details Checks if the stored pointers are equivalent. If either of the
   *  Elements is an empty pointer or the undefined element, result will always
   *  be false.
   *  \param l, r elements to compare.
   *  \return if the elements are the name. */
  inline bool operator==(Element l, Element r) {
    return (l && *l && r && *r) ? l == r : false;
  }
  
  // Inverse eq operators
  /*! \brief Equality test of integer and Element.
   *  \param l, r integer and Element to compare.
   *  \return if the element has the given atomic number.
   *  \see operator==(Element, uint8_t) */
  inline bool operator==(uint8_t l, Element r) { return r == l; }
  /*! \brief Equality test of strinf and Element.
   *  \param l, r string and Element to compare.
   *  \return if the element has the given atomic symbol or name.
   *  \see operator==(Element, std::string) */
  inline bool operator==(std::string l, Element r) { return r == l; }
  
  // Neq operators
  /*! \brief Inequality test of Element and integer.
   *  \param l, r Element and integer to compare.
   *  \return negation of equality test of the parameters.
   *  \see operator==(Element, uint8_t) */
  inline bool operator!=(Element l, uint8_t r) { return !(l == r); }
  /*! \brief Inequality test of integer and Element.
   *  \param l, r integer and Element to compare.
   *  \return negation of equality test of the parameters.
   *  \see operator==(Element, uint8_t) */
  inline bool operator!=(uint8_t l, Element r) { return !(r == l); }
  /*! \brief Inequality test of Element and string.
   *  \param l, r Element and string to compare.
   *  \return negation of equality test of the parameters.
   *  \see operator==(Element, std::string) */
  inline bool operator!=(Element l, std::string r) { return !(l == r); }
  /*! \brief Inequality test of string and Element.
   *  \param l, r string and Element to compare.
   *  \return negation of equality test of the parameters.
   *  \see operator==(Element, std::string) */
  inline bool operator!=(std::string l, Element r) { return !(r == l); }
  /*! \brief Inequality test of two Elements.
   *  \param l, r Elements to compare.
   *  \return negation of equality test of the parameters.
   *  \see operator==(Element, Element) */
  inline bool operator!=(Element l, Element r) { return !(l == r); }

} // namespace indigox

#endif /* INDIGOX_CLASSES_PERIODIC_TABLE_HPP */
