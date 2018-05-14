#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE indigoX
#include <boost/test/unit_test.hpp>

#include <indigox/utils/counter.hpp>
#include <indigox/utils/common.hpp>

using namespace indigox::utils;


BOOST_AUTO_TEST_CASE(ixcountableobject) {
  auto int_a = IXCountableObject<int>();
  auto int_b = IXCountableObject<int>();
  auto float_a = IXCountableObject<float>();
  auto float_b = IXCountableObject<float>();
  
  BOOST_TEST(int_a.GetUniqueID() + 1 == int_b.GetUniqueID());
  BOOST_TEST(float_a.GetUniqueID() + 1 == float_b.GetUniqueID());
  BOOST_TEST(int_a.GetUniqueID() == float_a.GetUniqueID());
}

