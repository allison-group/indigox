#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/common.hpp>
#include <indigox/utils/numerics.hpp>
#include <indigox/utils/doctest_proxy.hpp>

#ifndef INDIGOX_DISABLE_TESTS
#include <indigox/test/periodictable_test.hpp>
#endif

namespace indigox {
  
  test_suite_open("IXElement");
  
  IXElement::IXElement(uchar_ Z, string_ name, string_ sym, float_ mass,
                       uchar_ grp, uchar_ prd, uchar_ val, uchar_ oct,
                       uchar_ hyp, float_ rad, float_ cov, float_ vdw,
                       float_ chi)
  : _nme(name), _sym(sym), _grp(grp), _prd(prd), _Z(Z), _val(val), _oct(oct),
  _hyp(hyp), _mass(mass), _rad(rad), _cov(cov), _vdw(vdw), _chi(chi) { }
  
  test_case("IXElement construction") {
    test::TestElement e = test::CreateGenericTestElement();
    check_eq(e.get_Z(), 23);
    check_eq(e.get_nme(), "Testium");
    check_eq(e.get_sym(), "Tm");
    check_eq(e.get_mass(), approximately(12.34));
    check_eq(e.get_grp(), 5);
    check_eq(e.get_prd(), 6);
    check_eq(e.get_val(), 7);
    check_eq(e.get_oct(), 8);
    check_eq(e.get_hyp(), 9);
    check_eq(e.get_rad(), approximately(10.11));
    check_eq(e.get_cov(), approximately(12.13));
    check_eq(e.get_vdw(), approximately(14.1789));
    check_eq(e.get_chi(), approximately(15.8763));
  }
  
  string_ IXElement::ToString() const {
    std::stringstream ss;
    ss << _nme << "(" << +_Z << ", " << _sym << ")";
    return ss.str();
  }
  
  test_case("IXElement printing methods") {
    test::TestElement e = test::CreateGenericTestElement();
    check_eq(e.ToString(), "Testium(23, Tm)");
    std::stringstream ss;
    ss.str(""); ss << e.imp;
    check_eq(ss.str(), "Element(Testium)");
    ss.str(""); ss << *e.imp;
    check_eq(ss.str(), "Element(Testium)");
    ss.str(""); ss << Element();
    check_eq(ss.str(), "");
  }
  
  test_case("IXElement comparison methods") {
    Element e1 = test::CreateGenericTestElement().imp;
    Element e2 = test::TestElement(6, "Car", "C", 0.0, 5, 6, 7, 8, 9, 0.1,
                                   2.3, 4.9, 5.3).imp;
    
    // Check both shared_ptr comparisons
    check_ne(e1, e2);
    check_eq(e1, e1);
    // Check one direct, one shared_ptr comparison
    check_eq(e1, *e1);
    check_eq(*e1, e1);
    check_ne(e2, *e1);
    check_ne(*e2, e1);
    // Check both direct comparisons
    check_ne(*e1, *e2);
    check_eq(*e1, e1);
    // Check shared_ptr comparisions with integers (atomic number)
    check_eq(e1, 23);
    check_eq(23, e1);
    check_ne(23, e2);
    check_ne(e2, 23);
    // Check direct comparisions with integers (atomic number)
    check_eq(*e2, 6);
    check_eq(6, *e2);
    check_ne(6, *e1);
    check_ne(*e1, 6);
    // Check shared_ptr comparisons with strings (name)
    check_eq(e1, "Testium");
    check_eq(e1, "tesTIUM"); // name compare should be case insensitive
    check_eq("Testium", e1);
    check_eq("CAR", e2);     // name compare should be case insensitive
    check_ne("FAIL", e1);
    check_ne(e1, "FAIL");
    // Check shared_ptr comparisons with strings (symbol)
    check_eq(e1, "Tm");
    check_ne(e2, "Tm");
    check_ne(e1, "tm");  // symbol compare shoud be case sensitive
    check_eq("C", e2);
    check_ne("c", e2);   // symbol compare should be case sensitive
    // Check direct comparisons with strings (name)
    check_eq(*e1, "Testium");
    check_eq(*e1, "tesTIUM"); // name compare should be case insensitive
    check_eq("Testium", *e1);
    check_eq("CAR", *e2);
    check_ne("FAIL", *e1);
    check_ne(*e1, "FAIL");
    // Check direct comparisons with strings (symbol)
    check_eq(*e1, "Tm");
    check_ne(*e2, "Tm");
    check_ne(*e1, "tm");  // symbol compare shoud be case sensitive
    check_eq("C", *e2);
    check_ne("c", *e2);
  }
  
  test_case("IXElement property getting") {
    test::TestElement e = test::CreateGenericTestElement();
    check_eq(e.GetAtomicNumber(), 23);
    check_eq(e.GetName(), "Testium");
    check_eq(e.GetSymbol(), "Tm");
    check_eq(e.GetAtomicMass(), approximately(12.34));
    check_eq(e.GetGroup(), 5);
    check_eq(e.GetPeriod(), 6);
    check_eq(e.GetValenceElectronCount(), 7);
    check_eq(e.GetOctet(), 8);
    check_eq(e.GetHypervalentOctet(), 9);
    check_eq(e.GetAtomicRadius(), approximately(10.11));
    check_eq(e.GetCovalentRadius(), approximately(12.13));
    check_eq(e.GetVanDerWaalsRadius(), approximately(14.1789));
    check_eq(e.GetElectronegativity(), approximately(15.8763));
  }
  
  test_suite_close(); // IXElement
  
  //! \cond
  using erow = std::pair<Element, int>;
  
  std::ostream& operator<<(std::ostream& ss, erow er) {
    switch (er.second) {
      case 0:
        ss << "--- "; break;
      case 1:
        ss << std::setw(3) << er.first->GetAtomicNumber() << '|'; break;
      case 2:
        ss << std::setw(3) << er.first->GetSymbol() << '|'; break;
      case 3:
        ss << "    "; break;
      case 4:
      default:
        ss << "   |"; break;
    }
    return ss;
  }
  //! \endcond
  
  test_suite_open("IXPeriodicTable");
  
  string_ IXPeriodicTable::ToString() const {
    std::stringstream ss;
    size_ row_count = 0, restart = 0;
    std::vector<int_> elems = {
       1, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,-1,
       3, 4, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  6,  7,  8,  9, 10,-1,
      11,12, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 13, 14, 15, 16, 17, 18,-1,
      19,20,21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,-1,
      37,38,39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,-1,
      55,56,57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,-1,
      87,88,89,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,-1,
       0, 0, 0, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,  0,-1,
       0, 0, 0, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,  0,-1
    };
    if (!NumElements()) return "Empty PeriodicTable";
    ss << ' ';
    for (size_ i = 0; i < elems.size(); ) {
      if (elems[i] == -1) {
        ss << "\n";
        if ((row_count == 0 || row_count == 1) && elems[i-18] != 0) ss << '|';
        else ss << ' '; ++row_count;
        if (row_count < 3) i = restart;
        else { ++i; restart = i; row_count = 0; }
      } else if (elems[i] == 0 && i > 19 && elems[i-19] != 0 && row_count == 0) {
        ss << erow(_null, 0); ++i;
      } else if (elems[i] == 0 && elems[i+1] == 0) {
        ss << erow(_null, 3); ++i;
      } else if (elems[i] == 0 && elems[i+1] > 0 && row_count != 0) {
        ss << erow(_null, 4); ++i;
      } else if (elems[i] == 0 && elems[i+1] > 0 && row_count == 0) {
        ss << erow(_null, 3); ++i;
      } else if (elems[i] == 0 && elems[i+1] == -1) {
        ss << erow(_null, 3); ++i;
      } else {
        auto pos = _z_to.find(elems[i]);
        ss << erow(pos->second, row_count);
        ++i;
      }
    }
    // Final line at bottom of table
    for (size_ i = 0; i < 17; ++i) {
      if (i < 3) ss << erow(_null, 3);
      else ss << erow(_null, 0);
    }
    return ss.str();
  }
  
  test_case("IXPeriodicTable printing methods") {
    std::stringstream os;
    os << " ---                                                                 --- \n";
    os << "|  1|                                                               |  2|\n";
    os << "|  H|                                                               | He|\n";
    os << " --- ---                                         --- --- --- --- --- --- \n";
    os << "|  3|  4|                                       |  5|  6|  7|  8|  9| 10|\n";
    os << "| Li| Be|                                       |  B|  C|  N|  O|  F| Ne|\n";
    os << " --- ---                                         --- --- --- --- --- --- \n";
    os << "| 11| 12|                                       | 13| 14| 15| 16| 17| 18|\n";
    os << "| Na| Mg|                                       | Al| Si|  P|  S| Cl| Ar|\n";
    os << " --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n";
    os << "| 19| 20| 21| 22| 23| 24| 25| 26| 27| 28| 29| 30| 31| 32| 33| 34| 35| 36|\n";
    os << "|  K| Ca| Sc| Ti|  V| Cr| Mn| Fe| Co| Ni| Cu| Zn| Ga| Ge| As| Se| Br| Kr|\n";
    os << " --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n";
    os << "| 37| 38| 39| 40| 41| 42| 43| 44| 45| 46| 47| 48| 49| 50| 51| 52| 53| 54|\n";
    os << "| Rb| Sr|  Y| Zr| Nb| Mo| Tc| Ru| Rh| Pd| Ag| Cd| In| Sn| Sb| Te|  I| Xe|\n";
    os << " --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n";
    os << "| 55| 56| 57| 72| 73| 74| 75| 76| 77| 78| 79| 80| 81| 82| 83| 84| 85| 86|\n";
    os << "| Cs| Ba| La| Hf| Ta|  W| Re| Os| Ir| Pt| Au| Hg| Tl| Pb| Bi| Po| At| Rn|\n";
    os << " --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n";
    os << "| 87| 88| 89|104|105|106|107|108|109|110|111|112|113|114|115|116|117|118|\n";
    os << "| Fr| Ra| Ac| Db| Jl| Rf| Bh| Hn| Mt| Ds| Rg| Cn| Nh| Fl| Mc| Lv| Ts| Og|\n";
    os << " --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n";
    os << "            | 58| 59| 60| 61| 62| 63| 64| 65| 66| 67| 68| 69| 70| 71|    \n";
    os << "            | Ce| Pr| Nd| Pm| Sm| Eu| Gd| Tb| Dy| Ho| Er| Tm| Yb| Lu|    \n";
    os << "             --- --- --- --- --- --- --- --- --- --- --- --- --- ---     \n";
    os << "            | 90| 91| 92| 93| 94| 95| 96| 97| 98| 99|100|101|102|103|    \n";
    os << "            | Th| Pa|  U| Np| Pu| Am| Cm| Bk| Cf| Es| Fm| Md| No| Lr|    \n";
    os << "             --- --- --- --- --- --- --- --- --- --- --- --- --- --- ";
    string_ expected_to_string = os.str();
    
    test::TestPeriodicTable pt;
    os.str(""); os << PeriodicTable();
    check_eq(os.str(), "");
    os.str(""); os << pt.imp;
    check_eq(os.str(), "PeriodicTable(0 elements)");
    os.str(""); os << *pt.imp;
    check_eq(os.str(), "PeriodicTable(0 elements)");
    check_eq("Empty PeriodicTable", pt.ToString());
    
    pt.GeneratePeriodicTable();
    os.str(""); os << pt.imp;
    check_eq(os.str(), "PeriodicTable(118 elements)");
    os.str(""); os << *pt.imp;
    check_eq(os.str(), "PeriodicTable(118 elements)");
    check_eq(expected_to_string, pt.ToString());
  }
  
  Element IXPeriodicTable::GetElement(const uchar_ z) const noexcept {
    if (_z_to.find(z) != _z_to.end())
      return _z_to.at(z);
    return _null;
  }
  
  Element IXPeriodicTable::GetElement(const string_ name) const noexcept {
    std::string u;
    if (name.size() > 2) u = utils::ToUpperFirst(name);
    else u = name;
    if (_name_to.find(u) != _name_to.end()) {
      return _name_to.at(u);
    }
    return _null;
  }
  
  test_case("IXPeriodicTable getting elements") {
    // We build our test elements to compare the gets with
    test::TestPeriodicTable PT; PT.GeneratePeriodicTable();
    Element C = test::TestElement(6, "Carbon" , "C", 12.011000, 14, 2, 4, 8, 8,
                                  0.77, 0.77, 1.85, 2.55).imp;
    Element K = test::TestElement(19, "Potassium" , "K", 39.098300, 1, 4, 1, 2,
                                  2, 2.27, 2.03, 2.31, 0.82).imp;
    Element null_element = PT.get_null();
    
    check_eq(null_element, PT.GetUndefined());
    // Correct retrival
    check_eq(C, PT.GetElement(6));
    check_eq(K, PT.GetElement(19));
    check_eq(C, PT.GetElement("C"));
    check_eq(K, PT.GetElement("K"));
    check_eq(C, PT.GetElement("Carbon"));
    check_eq(K, PT.GetElement("Potassium"));
    // Symbol get case sensitive, name case insensitive
    check_ne(C, PT.GetElement("c"));
    check_eq(K, PT.GetElement("potASSIUM"));
    // Bad requests return null atom
    check_eq(null_element, PT.GetElement("k"));
    check_eq(null_element, PT.GetElement("NotAnElement"));
    check_eq(null_element, PT.GetElement(0));
    check_eq(null_element, PT.GetElement(255));
    // operators work
    check_eq(C, PT["C"]);
    check_eq(K, PT[19]);
    check_eq(C, PT["Carbon"]);
    check_eq(null_element, PT[255]);
    check_eq(null_element, PT["NotAnElement"]);
    
  }
  
#define INDIGOX_NUM_ELEMENTS 118
  void IXPeriodicTable::GeneratePeriodicTable() {
    _null = Element(   new IXElement(  0, "Undefined"    , "XX",   0.000000,  0, 0,  0, 0,  0,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(  1, new IXElement(  1, "Hydrogen"     , "H" ,   1.007970,  1, 1,  1, 2,  2,  0.78,  0.30,  1.20,  2.20));
    _z_to.emplace(  2, new IXElement(  2, "Helium"       , "He",   4.002602, 18, 1,  2, 2,  2,  1.28,  0.00,  1.22,  0.00));
    _z_to.emplace(  3, new IXElement(  3, "Lithium"      , "Li",   6.941000,  1, 2,  1, 2,  2,  1.52,  1.23,  0.00,  0.98));
    _z_to.emplace(  4, new IXElement(  4, "Beryllium"    , "Be",   9.012182,  2, 2,  2, 4,  4,  1.13,  0.89,  0.00,  1.57));
    _z_to.emplace(  5, new IXElement(  5, "Boron"        , "B" ,  10.811000, 13, 2,  3, 6,  6,  0.83,  0.88,  2.08,  2.04));
    _z_to.emplace(  6, new IXElement(  6, "Carbon"       , "C" ,  12.011000, 14, 2,  4, 8,  8,  0.77,  0.77,  1.85,  2.55));
    _z_to.emplace(  7, new IXElement(  7, "Nitrogen"     , "N" ,  14.006740, 15, 2,  5, 8,  8,  0.71,  0.70,  1.54,  3.04));
    _z_to.emplace(  8, new IXElement(  8, "Oxygen"       , "O" ,  15.999400, 16, 2,  6, 8,  8,  0.60,  0.66,  1.40,  3.44));
    _z_to.emplace(  9, new IXElement(  9, "Fluorine"     , "F" ,  18.998403, 17, 2,  7, 8,  8,  0.71,  0.58,  1.35,  3.98));
    _z_to.emplace( 10, new IXElement( 10, "Neon"         , "Ne",  20.179700, 18, 2,  8, 8,  8,  0.00,  0.00,  1.60,  0.00));
    _z_to.emplace( 11, new IXElement( 11, "Sodium"       , "Na",  22.989768,  1, 3,  1, 2,  2,  1.54,  0.00,  2.31,  0.93));
    _z_to.emplace( 12, new IXElement( 12, "Magnesium"    , "Mg",  24.305060,  2, 3,  2, 4,  4,  1.60,  1.36,  0.00,  1.31));
    _z_to.emplace( 13, new IXElement( 13, "Aluminum"     , "Al",  26.981539, 13, 3,  3, 6,  6,  1.43,  1.25,  2.05,  1.61));
    _z_to.emplace( 14, new IXElement( 14, "Silicon"      , "Si",  28.085500, 14, 3,  4, 8,  8,  1.17,  1.17,  2.00,  1.90));
    _z_to.emplace( 15, new IXElement( 15, "Phosphorus"   , "P" ,  30.973762, 15, 3,  5, 8, 10,  1.15,  1.10,  1.90,  2.19));
    _z_to.emplace( 16, new IXElement( 16, "Sulfur"       , "S" ,  32.066000, 16, 3,  6, 8, 12,  1.04,  1.04,  1.85,  2.58));
    _z_to.emplace( 17, new IXElement( 17, "Chlorine"     , "Cl",  35.452700, 17, 3,  7, 8,  8,  0.00,  0.99,  1.81,  3.16));
    _z_to.emplace( 18, new IXElement( 18, "Argon"        , "Ar",  39.948000, 18, 3,  8, 8,  8,  1.74,  0.00,  1.91,  0.00));
    _z_to.emplace( 19, new IXElement( 19, "Potassium"    , "K" ,  39.098300,  1, 4,  1, 2,  2,  2.27,  2.03,  2.31,  0.82));
    _z_to.emplace( 20, new IXElement( 20, "Calcium"      , "Ca",  40.078000,  2, 4,  2, 4,  4,  1.97,  1.74,  0.00,  1.00));
    _z_to.emplace( 21, new IXElement( 21, "Scandium"     , "Sc",  44.955910,  3, 4,  3, 8,  8,  1.61,  1.44,  0.00,  1.36));
    _z_to.emplace( 22, new IXElement( 22, "Titanium"     , "Ti",  47.867000,  4, 4,  4, 8,  8,  1.45,  1.32,  0.00,  1.54));
    _z_to.emplace( 23, new IXElement( 23, "Vanadium"     , "V" ,  50.941500,  5, 4,  5, 8,  8,  1.32,  0.00,  0.00,  1.63));
    _z_to.emplace( 24, new IXElement( 24, "Chromium"     , "Cr",  51.996100,  6, 4,  6, 8,  8,  1.25,  0.00,  0.00,  1.66));
    _z_to.emplace( 25, new IXElement( 25, "Manganese"    , "Mn",  54.938050,  7, 4,  7, 8,  8,  1.24,  1.77,  0.00,  1.55));
    _z_to.emplace( 26, new IXElement( 26, "Iron"         , "Fe",  55.845000,  8, 4,  8, 8,  8,  1.24,  1.16,  0.00,  1.83));
    _z_to.emplace( 27, new IXElement( 27, "Cobalt"       , "Co",  58.933200,  9, 4,  9, 8,  8,  1.25,  1.16,  0.00,  1.88));
    _z_to.emplace( 28, new IXElement( 28, "Nickel"       , "Ni",  58.693400, 10, 4, 10, 8,  8,  1.25,  1.15,  0.00,  1.91));
    _z_to.emplace( 29, new IXElement( 29, "Copper"       , "Cu",  63.546000, 11, 4, 11, 8,  8,  1.28,  1.17,  0.00,  1.90));
    _z_to.emplace( 30, new IXElement( 30, "Zinc"         , "Zn",  65.390000, 12, 4, 12, 8,  8,  1.33,  1.25,  0.00,  1.65));
    _z_to.emplace( 31, new IXElement( 31, "Gallium"      , "Ga",  69.723000, 13, 4,  3, 6,  6,  1.22,  1.25,  0.00,  1.81));
    _z_to.emplace( 32, new IXElement( 32, "Germanium"    , "Ge",  72.610000, 14, 4,  4, 8,  8,  1.23,  1.22,  0.00,  2.01));
    _z_to.emplace( 33, new IXElement( 33, "Arsenic"      , "As",  74.921590, 15, 4,  5, 8,  8,  1.25,  1.21,  2.00,  2.18));
    _z_to.emplace( 34, new IXElement( 34, "Selenium"     , "Se",  78.960000, 16, 4,  6, 8,  8,  2.15,  1.17,  2.00,  2.55));
    _z_to.emplace( 35, new IXElement( 35, "Bromine"      , "Br",  79.904000, 17, 4,  7, 8, 12,  0.00,  1.14,  1.95,  2.96));
    _z_to.emplace( 36, new IXElement( 36, "Krypton"      , "Kr",  83.800000, 18, 4,  8, 8,  8,  0.00,  1.89,  1.98,  0.00));
    _z_to.emplace( 37, new IXElement( 37, "Rubidium"     , "Rb",  85.467800,  1, 5,  1, 2,  2,  1.48,  0.00,  2.44,  0.82));
    _z_to.emplace( 38, new IXElement( 38, "Strontium"    , "Sr",  87.620000,  2, 5,  2, 4,  4,  2.15,  1.92,  0.00,  0.95));
    _z_to.emplace( 39, new IXElement( 39, "Yttrium"      , "Y" ,  88.905850,  3, 5,  3, 8,  8,  1.81,  1.62,  0.00,  1.22));
    _z_to.emplace( 40, new IXElement( 40, "Zirconium"    , "Zr",  91.224000,  4, 5,  4, 8,  8,  1.60,  1.45,  0.00,  1.30));
    _z_to.emplace( 41, new IXElement( 41, "Niobium"      , "Nb",  92.906380,  5, 5,  5, 8,  8,  1.43,  1.34,  0.00,  1.60));
    _z_to.emplace( 42, new IXElement( 42, "Molybdenum"   , "Mo",  95.940000,  6, 5,  6, 8,  8,  1.36,  1.29,  0.00,  2.16));
    _z_to.emplace( 43, new IXElement( 43, "Technetium"   , "Tc",  98.907200,  7, 5,  7, 8,  8,  1.36,  0.00,  0.00,  1.90));
    _z_to.emplace( 44, new IXElement( 44, "Ruthenium"    , "Ru", 101.070000,  8, 5,  8, 8,  8,  1.34,  1.24,  0.00,  2.20));
    _z_to.emplace( 45, new IXElement( 45, "Rhodium"      , "Rh", 102.905500,  9, 5,  9, 8,  8,  1.34,  1.25,  0.00,  2.28));
    _z_to.emplace( 46, new IXElement( 46, "Palladium"    , "Pd", 106.420000, 10, 5, 10, 8,  8,  1.38,  1.28,  0.00,  2.20));
    _z_to.emplace( 47, new IXElement( 47, "Silver"       , "Ag", 107.868200, 11, 5, 11, 8,  8,  1.44,  1.34,  0.00,  1.93));
    _z_to.emplace( 48, new IXElement( 48, "Cadmium"      , "Cd", 112.411000, 12, 5, 12, 8,  8,  1.49,  1.41,  0.00,  1.69));
    _z_to.emplace( 49, new IXElement( 49, "Indium"       , "In", 114.818000, 13, 5,  3, 6,  6,  1.63,  1.50,  0.00,  1.78));
    _z_to.emplace( 50, new IXElement( 50, "Tin"          , "Sn", 118.710000, 14, 5,  4, 8,  8,  1.41,  1.40,  2.00,  1.96));
    _z_to.emplace( 51, new IXElement( 51, "Antimony"     , "Sb", 121.760000, 15, 5,  5, 8,  8,  1.82,  1.41,  2.20,  2.05));
    _z_to.emplace( 52, new IXElement( 52, "Tellurium"    , "Te", 127.600000, 16, 5,  6, 8,  8,  1.43,  1.37,  2.20,  2.10));
    _z_to.emplace( 53, new IXElement( 53, "Iodine"       , "I" , 126.904470, 17, 5,  7, 8,  8,  0.00,  1.33,  2.15,  2.66));
    _z_to.emplace( 54, new IXElement( 54, "Xenon"        , "Xe", 131.290000, 18, 5,  8, 8,  8,  2.18,  2.09,  2.16,  2.60));
    _z_to.emplace( 55, new IXElement( 55, "Caesium"      , "Cs", 132.905430,  1, 6,  1, 2,  2,  2.65,  2.35,  2.62,  0.79));
    _z_to.emplace( 56, new IXElement( 56, "Barium"       , "Ba", 137.327000,  2, 6,  2, 4,  4,  2.17,  1.98,  0.00,  0.89));
    _z_to.emplace( 57, new IXElement( 57, "Lanthanum"    , "La", 138.905500,  3, 6,  3, 8,  8,  1.88,  1.69,  0.00,  1.10));
    _z_to.emplace( 58, new IXElement( 58, "Cerium"       , "Ce", 140.115000,  0, 6,  0, 8,  8,  1.82,  1.65,  0.00,  1.12));
    _z_to.emplace( 59, new IXElement( 59, "Praseodymium" , "Pr", 140.907650,  0, 6,  0, 8,  8,  1.83,  1.65,  0.00,  1.13));
    _z_to.emplace( 60, new IXElement( 60, "Neodymium"    , "Nd", 144.240000,  0, 6,  0, 8,  8,  1.82,  1.64,  0.00,  1.14));
    _z_to.emplace( 61, new IXElement( 61, "Promethium"   , "Pm", 144.912700,  0, 6,  0, 8,  8,  1.81,  0.00,  0.00,  0.94));
    _z_to.emplace( 62, new IXElement( 62, "Samarium"     , "Sm", 150.360000,  0, 6,  0, 8,  8,  1.80,  1.66,  0.00,  1.17));
    _z_to.emplace( 63, new IXElement( 63, "Europium"     , "Eu", 151.965000,  0, 6,  0, 8,  8,  2.04,  1.85,  0.00,  1.20));
    _z_to.emplace( 64, new IXElement( 64, "Gadolinium"   , "Gd", 157.250000,  0, 6,  0, 8,  8,  1.80,  1.61,  0.00,  0.94));
    _z_to.emplace( 65, new IXElement( 65, "Terbium"      , "Tb", 158.925340,  0, 6,  0, 8,  8,  1.78,  1.59,  0.00,  1.22));
    _z_to.emplace( 66, new IXElement( 66, "Dysprosium"   , "Dy", 162.500000,  0, 6,  0, 8,  8,  1.77,  1.59,  0.00,  1.23));
    _z_to.emplace( 67, new IXElement( 67, "Holmium"      , "Ho", 164.930320,  0, 6,  0, 8,  8,  1.77,  1.58,  0.00,  1.24));
    _z_to.emplace( 68, new IXElement( 68, "Erbium"       , "Er", 167.260000,  0, 6,  0, 8,  8,  1.76,  1.57,  0.00,  1.25));
    _z_to.emplace( 69, new IXElement( 69, "Thulium"      , "Tm", 168.934210,  0, 6,  0, 8,  8,  1.75,  1.56,  0.00,  0.96));
    _z_to.emplace( 70, new IXElement( 70, "Ytterbium"    , "Yb", 173.040000,  0, 6,  0, 8,  8,  1.94,  1.70,  0.00,  1.27));
    _z_to.emplace( 71, new IXElement( 71, "Lutetium"     , "Lu", 174.967000,  3, 6,  3, 8,  8,  1.72,  1.56,  0.00,  1.30));
    _z_to.emplace( 72, new IXElement( 72, "Hafnium"      , "Hf", 178.490000,  4, 6,  4, 8,  8,  1.56,  1.44,  0.00,  1.50));
    _z_to.emplace( 73, new IXElement( 73, "Tantalum"     , "Ta", 180.947900,  5, 6,  5, 8,  8,  1.43,  1.34,  0.00,  2.36));
    _z_to.emplace( 74, new IXElement( 74, "Tungsten"     , "W" , 183.840000,  6, 6,  6, 8,  8,  1.37,  1.30,  0.00,  1.90));
    _z_to.emplace( 75, new IXElement( 75, "Rhenium"      , "Re", 186.207000,  7, 6,  7, 8,  8,  1.37,  1.28,  0.00,  2.20));
    _z_to.emplace( 76, new IXElement( 76, "Osmium"       , "Os", 190.230000,  8, 6,  8, 8,  8,  1.35,  1.26,  0.00,  2.20));
    _z_to.emplace( 77, new IXElement( 77, "Iridium"      , "Ir", 192.217000,  9, 6,  9, 8,  8,  1.36,  1.26,  0.00,  2.28));
    _z_to.emplace( 78, new IXElement( 78, "Platinum"     , "Pt", 195.080000, 10, 6, 10, 8,  8,  1.38,  1.29,  0.00,  2.54));
    _z_to.emplace( 79, new IXElement( 79, "Gold"         , "Au", 196.966540, 11, 6, 11, 8,  8,  1.44,  1.34,  0.00,  2.00));
    _z_to.emplace( 80, new IXElement( 80, "Mercury"      , "Hg", 200.590000, 12, 6, 12, 8,  8,  1.60,  1.44,  0.00,  1.80));
    _z_to.emplace( 81, new IXElement( 81, "Thallium"     , "Tl", 204.383300, 13, 6,  3, 6,  6,  1.70,  1.55,  0.00,  2.33));
    _z_to.emplace( 82, new IXElement( 82, "Lead"         , "Pb", 207.200000, 14, 6,  4, 8,  8,  1.75,  1.54,  0.00,  2.02));
    _z_to.emplace( 83, new IXElement( 83, "Bismuth"      , "Bi", 208.980370, 15, 6,  5, 8,  8,  1.55,  1.52,  2.40,  2.00));
    _z_to.emplace( 84, new IXElement( 84, "Polonium"     , "Po", 208.982400, 16, 6,  6, 8,  8,  1.67,  1.53,  0.00,  2.20));
    _z_to.emplace( 85, new IXElement( 85, "Astatine"     , "At", 209.987100, 17, 6,  7, 8,  8,  0.00,  0.00,  0.00,  1.96));
    _z_to.emplace( 86, new IXElement( 86, "Radon"        , "Rn", 222.017600, 18, 6,  8, 8,  8,  0.00,  0.00,  0.00,  0.70));
    _z_to.emplace( 87, new IXElement( 87, "Francium"     , "Fr", 223.019700,  1, 7,  1, 2,  2,  2.70,  0.00,  0.00,  0.70));
    _z_to.emplace( 88, new IXElement( 88, "Radium"       , "Ra", 226.025400,  2, 7,  2, 4,  4,  2.23,  0.00,  0.00,  0.89));
    _z_to.emplace( 89, new IXElement( 89, "Actinium"     , "Ac", 227.027800,  3, 7,  3, 8,  8,  1.88,  0.00,  0.00,  1.30));
    _z_to.emplace( 90, new IXElement( 90, "Thorium"      , "Th", 232.038100,  0, 7,  0, 8,  8,  1.80,  0.00,  0.00,  0.00));
    _z_to.emplace( 91, new IXElement( 91, "Protactinium" , "Pa", 231.035880,  0, 7,  0, 8,  8,  1.61,  0.00,  0.00,  1.38));
    _z_to.emplace( 92, new IXElement( 92, "Uranium"      , "U" , 238.028900,  0, 7,  0, 8,  8,  1.54,  0.00,  0.00,  1.26));
    _z_to.emplace( 93, new IXElement( 93, "Neptunium"    , "Np", 237.048200,  0, 7,  0, 8,  8,  1.50,  0.00,  0.00,  1.28));
    _z_to.emplace( 94, new IXElement( 94, "Plutonium"    , "Pu", 244.064200,  7, 0,  7, 8,  8,  0.00,  0.00,  0.00,  1.30));
    _z_to.emplace( 95, new IXElement( 95, "Americium"    , "Am", 243.061400,  0, 7,  0, 8,  8,  1.73,  0.00,  0.00,  1.30));
    _z_to.emplace( 96, new IXElement( 96, "Curium"       , "Cm", 247.070300,  0, 7,  0, 8,  8,  1.74,  0.00,  0.00,  1.30));
    _z_to.emplace( 97, new IXElement( 97, "Berkelium"    , "Bk", 247.070300,  0, 7,  0, 8,  8,  1.70,  0.00,  0.00,  1.30));
    _z_to.emplace( 98, new IXElement( 98, "Californium"  , "Cf", 251.079600,  0, 7,  0, 8,  8,  1.69,  0.00,  0.00,  1.30));
    _z_to.emplace( 99, new IXElement( 99, "Einsteinium"  , "Es", 252.083000,  0, 7,  0, 8,  8,  2.03,  0.00,  0.00,  1.30));
    _z_to.emplace(100, new IXElement(100, "Fermium"      , "Fm", 257.095100,  0, 7,  0, 8,  8,  0.00,  0.00,  0.00,  1.30));
    _z_to.emplace(101, new IXElement(101, "Mendelevium"  , "Md", 258.100000,  0, 7,  0, 8,  8,  0.00,  0.00,  0.00,  1.30));
    _z_to.emplace(102, new IXElement(102, "Nobelium"     , "No", 259.100900,  0, 7,  0, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(103, new IXElement(103, "Lawrencium"   , "Lr", 262.110000,  3, 7,  3, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(104, new IXElement(104, "Dubnium"      , "Db", 261.110000,  4, 7,  4, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(105, new IXElement(105, "Joliotium"    , "Jl", 262.114000,  5, 7,  5, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(106, new IXElement(106, "Rutherfordium", "Rf", 263.118000,  6, 7,  6, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(107, new IXElement(107, "Bohrium"      , "Bh", 262.120000,  7, 7,  7, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(108, new IXElement(108, "Hahnium"      , "Hn",   0.000000,  8, 7,  8, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(109, new IXElement(109, "Meitnerium"   , "Mt",   0.000000,  9, 7,  9, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(110, new IXElement(110, "Darmstadtium" , "Ds",   0.000000, 10, 7, 10, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(111, new IXElement(111, "Roentgenium"  , "Rg",   0.000000, 11, 7, 11, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(112, new IXElement(112, "Copernicium"  , "Cn",   0.000000, 12, 7, 12, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(113, new IXElement(113, "Nihonium"     , "Nh",   0.000000, 13, 7,  3, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(114, new IXElement(114, "Flerovium"    , "Fl",   0.000000, 14, 7,  4, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(115, new IXElement(115, "Moscovium"    , "Mc",   0.000000, 15, 7,  5, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(116, new IXElement(116, "Livermorium"  , "Lv",   0.000000, 16, 7,  6, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(117, new IXElement(117, "Tennessine"   , "Ts",   0.000000, 17, 7,  7, 8,  8,  0.00,  0.00,  0.00,  0.00));
    _z_to.emplace(118, new IXElement(118, "Oganesson"    , "Og",   0.000000, 18, 7,  8, 8,  8,  0.00,  0.00,  0.00,  0.00));
    
    for (auto& ze: _z_to) {
      _name_to.emplace(ze.second->GetSymbol(), ze.second);
      _name_to.emplace(ze.second->GetName(), ze.second);
    }
  }
  
  test_case("IXPeriodicTable generation") {
    test::TestPeriodicTable PT; PT.GeneratePeriodicTable();
    test::TestPeriodicTable::ZType& z_to_element = PT.get_z_to();
    test::TestPeriodicTable::XType& name_to_element = PT.get_name_to();
    
    string_ name;
    uchar_ z;
    Element e;
    subcase("Checking atomic number table") {
      check_eq(z_to_element.size(), INDIGOX_NUM_ELEMENTS);
      DOCTEST_INFO("Current atomic number:");
      for (auto& ze : z_to_element) {
        std::tie(z, e) = ze;
        DOCTEST_CAPTURE(z);
        check_eq(e, name_to_element.at(e->GetSymbol()));
        check_eq(e, name_to_element.at(e->GetName()));
      }
    }
    
    subcase("Checking atomic name/symbol table") {
      check_eq(name_to_element.size(), INDIGOX_NUM_ELEMENTS * 2);
      DOCTEST_INFO("Current atomic symbol/name:");
      for (auto& ne : name_to_element) {
        std::tie(name, e) = ne;
        DOCTEST_CAPTURE(name);
        check_eq(e, z_to_element.at(e->GetAtomicNumber()));
      }
    }
    
    subcase("Checking null atom") {
      IXElement null_ = *PT.get_null();
      check_eq(null_.GetName(), "Undefined");
      check_eq(null_.GetSymbol(), "XX");
      check_eq(null_.GetAtomicNumber(), 0);
      check_eq(PT.NumElements(), INDIGOX_NUM_ELEMENTS);
    }
    
  }
  
  PeriodicTable GetPeriodicTable() {
    static PeriodicTable instance = PeriodicTable();
    if (!instance) {
      instance.reset(new IXPeriodicTable());
      instance->GeneratePeriodicTable();
    }
    return instance;
  }
  
  test_case("IXPeriodicTable GetPeriodicTable helper method") {
    PeriodicTable PT = GetPeriodicTable();
    check_eq(PT->NumElements(), INDIGOX_NUM_ELEMENTS);
    PeriodicTable PT2 = GetPeriodicTable();
    check_eq(PT2->NumElements(), INDIGOX_NUM_ELEMENTS);
    check_eq(PT.get(), PT2.get());
  }
  
  test_suite_close(); // IXPeriodicTable
  
} // namespace indigox

