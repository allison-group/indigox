#define BOOST_TEST_MODULE PeriodicTable test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/helpers.hpp>

#include <algorithm>
#include <iostream>

using namespace indigox;
// Test for the correct generation of PeriodicTable instances
BOOST_AUTO_TEST_CASE(periodictable_instance_test,
                     *boost::unit_test::expected_failures(1)) {
  PeriodicTable pt = GetPeriodicTable();
  PeriodicTable pt2 = GetPeriodicTable();
  
  BOOST_TEST(pt->NumElements() == 118);
  BOOST_TEST(pt == pt2);
  string_ expected_table = "\
 ---                                                                 --- \n\
|  1|                                                               |  2|\n\
|  H|                                                               | He|\n\
 --- ---                                         --- --- --- --- --- --- \n\
|  3|  4|                                       |  5|  6|  7|  8|  9| 10|\n\
| Li| Be|                                       |  B|  C|  N|  O|  F| Ne|\n\
 --- ---                                         --- --- --- --- --- --- \n\
| 11| 12|                                       | 13| 14| 15| 16| 17| 18|\n\
| Na| Mg|                                       | Al| Si|  P|  S| Cl| Ar|\n\
 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n\
| 19| 20| 21| 22| 23| 24| 25| 26| 27| 28| 29| 30| 31| 32| 33| 34| 35| 36|\n\
|  K| Ca| Sc| Ti|  V| Cr| Mn| Fe| Co| Ni| Cu| Zn| Ga| Ge| As| Se| Br| Kr|\n\
 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n\
| 37| 38| 39| 40| 41| 42| 43| 44| 45| 46| 47| 48| 49| 50| 51| 52| 53| 54|\n\
| Rb| Sr|  Y| Zr| Nb| Mo| Tc| Ru| Rh| Pd| Ag| Cd| In| Sn| Sb| Te|  I| Xe|\n\
 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n\
| 55| 56| 57| 72| 73| 74| 75| 76| 77| 78| 79| 80| 81| 82| 83| 84| 85| 86|\n\
| Cs| Ba| La| Hf| Ta|  W| Re| Os| Ir| Pt| Au| Hg| Tl| Pb| Bi| Po| At| Rn|\n\
 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n\
| 87| 88| 89|104|105|106|107|108|109|110|111|112|113|114|115|116|117|118|\n\
| Fr| Ra| Ac| Db| Jl| Rf| Bh| Hn| Mt| Ds| Rg| Cn| Nh| Fl| Mc| Lv| Ts| Og|\n\
 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n\
            | 58| 59| 60| 61| 62| 63| 64| 65| 66| 67| 68| 69| 70| 71|\n\
            | Ce| Pr| Nd| Pm| Sm| Eu| Gd| Tb| Dy| Ho| Er| Tm| Yb| Lu|\n\
             --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n\
            | 90| 91| 92| 93| 94| 95| 96| 97| 98| 99|100|101|102|103|\n\
            | Th| Pa|  U| Np| Pu| Am| Cm| Bk| Cf| Es| Fm| Md| No| Lr|\n\
             --- --- --- --- --- --- --- --- --- --- --- --- --- --- \n";
  
  BOOST_CHECK_MESSAGE(pt->ToString() == expected_table,
                      "Full periodic table string does not match expected.");
}

// Test for correct handling of the null element
BOOST_AUTO_TEST_CASE(null_element_test) {
  PeriodicTable pt = GetPeriodicTable();
  Element null_e = pt->GetUndefinedElement();
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

namespace bdata = boost::unit_test::data;

// Test for correct retrieval of elements
BOOST_DATA_TEST_CASE(element_retrieve_test,
                     bdata::random(1, 118) ^ bdata::xrange(10),
                     random_num, index) {
  // Test retrieving 10 random elements
  PeriodicTable pt = GetPeriodicTable();
  Element e1 = pt->GetElement(random_num);
  BOOST_TEST(e1);
  // Check correct element obtained
  BOOST_TEST(e1->GetAtomicNumber() == random_num);
  // Check other getting methods return same element
  BOOST_TEST(e1 == (*pt)[random_num]);
  BOOST_TEST(e1 == (*pt)[e1->GetSymbol()]);
  BOOST_TEST(e1 == (*pt)[e1->GetName()]);
  BOOST_TEST(e1 == pt->GetElement(e1->GetSymbol()));
  BOOST_TEST(e1 == pt->GetElement(e1->GetName()));
  // Check that the name and symbol are correctly size;
  BOOST_TEST(e1->GetSymbol().size() <= 2);
  BOOST_TEST(e1->GetName().size() > 2);
}

// Test comparison operators for elements
BOOST_AUTO_TEST_CASE(element_compare_test) {
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
