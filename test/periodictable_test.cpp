#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/helpers.hpp>

#include <algorithm>
#include <iostream>
#include <sstream>

BOOST_AUTO_TEST_SUITE(ixperiodictable);

using namespace indigox;
// Test for the correct generation of PeriodicTable instances
BOOST_AUTO_TEST_CASE(instance_constructor,
                     *boost::unit_test::expected_failures(1)) {
  PeriodicTable pt = GetPeriodicTable();
  PeriodicTable pt2 = GetPeriodicTable();
  
  BOOST_TEST(pt->NumElements() == 118);
  BOOST_TEST(pt == pt2);
}

// Test for correct handling of the null element
BOOST_AUTO_TEST_CASE(null_element) {
  PeriodicTable pt = GetPeriodicTable();
  Element null_e = pt->GetUndefined();
  // Check that the null element is not an empty pointer
  BOOST_TEST(null_e);
  BOOST_TEST(bool(*null_e) == false);
  // Check the symbol, num, name correctly defined
  BOOST_TEST(null_e->GetAtomicNumber() == 0);
  BOOST_TEST(null_e->GetSymbol() == "XX");
  BOOST_TEST(null_e->GetName() == "Undefined");
  // Check null_e always compares false
  BOOST_TEST(null_e != null_e);
  BOOST_TEST(null_e != 0);
  BOOST_TEST(null_e != "XX");
  BOOST_TEST(null_e != "Undefined");
}


// Test for correct retrieval of elements
BOOST_AUTO_TEST_CASE(element_retrieve) {
  PeriodicTable pt = GetPeriodicTable();
  Element carbon = pt->GetElement(6);
  BOOST_TEST((carbon && *carbon));
  // Check correct element obtained
  BOOST_TEST(carbon->GetAtomicNumber() == 6);
  // Check other getting methods return same element
  BOOST_TEST(carbon == (*pt)[6]);
  BOOST_TEST(carbon == (*pt)[carbon->GetSymbol()]);
  BOOST_TEST(carbon == (*pt)[carbon->GetName()]);
  BOOST_TEST(carbon == pt->GetElement(carbon->GetSymbol()));
  BOOST_TEST(carbon == pt->GetElement(carbon->GetName()));
  BOOST_TEST(carbon == pt->GetElement("caRBOn"));  // name get case insensitive
  // Check that the name and symbol are correct
  BOOST_TEST(carbon->GetSymbol() == "C");
  BOOST_TEST(carbon->GetName() == "Carbon");
}

// Test comparison operators for elements
BOOST_AUTO_TEST_CASE(element_compare) {
  PeriodicTable pt = GetPeriodicTable();
  Element carbon = pt->GetElement(6);
  Element nitrogen = pt->GetElement(7);
  BOOST_TEST(carbon == carbon);
  BOOST_TEST(carbon == 6);
  BOOST_TEST(6 == carbon);
  BOOST_TEST(carbon == "C");
  BOOST_TEST("C" == carbon);
  BOOST_TEST(carbon == "Carbon");
  BOOST_TEST("Carbon" == carbon);
  BOOST_TEST(carbon == "caRBoN"); // Name compare should be case insensitive
  BOOST_TEST("caRBoN" == carbon);
  BOOST_TEST(carbon != nitrogen);
  BOOST_TEST(nitrogen != carbon);
  BOOST_TEST(carbon != 7);
  BOOST_TEST(7 != carbon);
  BOOST_TEST(carbon != "N");
  BOOST_TEST("N" != carbon);
  BOOST_TEST(carbon != "c");  // Symbol compare should be case sensitive
  BOOST_TEST("c" != carbon);
  BOOST_TEST(carbon != "Nitrogen");
  BOOST_TEST("Nitrogen" != carbon);
}

// Test bad requests handling
BOOST_AUTO_TEST_CASE(bad_element_retrieve) {
  PeriodicTable pt = GetPeriodicTable();
  // Expect no throw at the extremes
  BOOST_CHECK_NO_THROW(pt->GetElement(1));
  BOOST_CHECK_NO_THROW(pt->GetElement(pt->NumElements()));
  // Expect invalid_argument when out of range
  BOOST_CHECK_THROW(pt->GetElement(0), std::invalid_argument);
  BOOST_CHECK_THROW(pt->GetElement(pt->NumElements()+1), std::invalid_argument);
  // Expect no throw when getting valid elements
  BOOST_CHECK_NO_THROW(pt->GetElement("C"));
  BOOST_CHECK_NO_THROW(pt->GetElement("Cl"));
  // Expect invalid_argument when getting invalid name/symbol
  BOOST_CHECK_THROW(pt->GetElement("Xx"), std::invalid_argument);
  BOOST_CHECK_THROW(pt->GetElement("NotARealElementium"), std::invalid_argument);
  // Expect case sensitive getting of symbol, but not name
  BOOST_CHECK_THROW(pt->GetElement("c"), std::invalid_argument);
  BOOST_CHECK_NO_THROW(pt->GetElement("carbon"));
}

// Test ostreams
BOOST_AUTO_TEST_CASE(printing_methods) {
  PeriodicTable pt = GetPeriodicTable();
  PeriodicTable pt_fail;
  boost::test_tools::output_test_stream os;
  os << pt;
  BOOST_TEST(os.is_equal("PeriodicTable(118 elements)"));
  os << pt_fail;
  BOOST_TEST(os.is_empty());
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
  
  BOOST_TEST(os.is_equal(pt->ToString()));
  
  os << Element();
  BOOST_TEST(os.is_empty());
  os << pt->GetElement(6);
  BOOST_TEST(os.is_equal("Element(Carbon)"));
  os << "6-Carbon (C)";
  BOOST_TEST(os.is_equal(pt->GetElement(6)->ToString()));
}

// Test element properties get
BOOST_AUTO_TEST_CASE(element_get_properties) {
  PeriodicTable pt = GetPeriodicTable();
  Element cm = pt->GetElement("Cm"), w = pt->GetElement("W");
  BOOST_TEST((cm->GetAtomicMass() == 247.0703 && w->GetAtomicMass() == 183.84),
             boost::test_tools::tolerance(0.0000001));
  BOOST_TEST((cm->GetAtomicNumber() == 96 && w->GetAtomicNumber() == 74));
  BOOST_TEST((cm->GetName() == "Curium" && w->GetName() == "Tungsten"));
  BOOST_TEST((cm->GetSymbol() == "Cm" && w->GetSymbol() == "W"));
  BOOST_TEST((cm->GetGroup() == 0 && w->GetGroup() == 6));
  BOOST_TEST((cm->GetPeriod() == 7 && w->GetPeriod() == 6));
  BOOST_TEST((cm->GetValenceElectronCount() == 0 && w->GetValenceElectronCount() == 6));
  BOOST_TEST((cm->GetOctet() == 8 && w->GetOctet() == 8));
  BOOST_TEST((cm->GetHypervalentOctet() == 8 && w->GetHypervalentOctet() == 8));
  BOOST_TEST((cm->GetAtomicRadius() == 1.74 && w->GetAtomicRadius() == 1.37),
             boost::test_tools::tolerance(0.0000001));
  BOOST_TEST((cm->GetCovalentRadius() == 0.00 && w->GetCovalentRadius() == 1.30),
             boost::test_tools::tolerance(0.0000001));
  BOOST_TEST((cm->GetVanDerWaalsRadius() == 0.00 && w->GetVanDerWaalsRadius() == 0.00),
             boost::test_tools::tolerance(0.0000001));
  BOOST_TEST((cm->GetElectronegativity() == 1.30 && w->GetElectronegativity() == 1.90),
             boost::test_tools::tolerance(0.0000001));
}
BOOST_AUTO_TEST_SUITE_END();
