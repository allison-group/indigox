#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <indigox/utils/doctest_proxy.hpp>

#include <indigox/utils/common.hpp>
#include <indigox/utils/counter.hpp>

//test_case("IXCountableObject counts correctly") {
//  auto int_a = indigox::utils::IXCountableObject<int>();
//  auto int_b = indigox::utils::IXCountableObject<int>();
//  auto float_a = indigox::utils::IXCountableObject<float>();
//  auto float_b = indigox::utils::IXCountableObject<float>();
//  
//  check(int_a.GetUniqueID() + 1 == int_b.GetUniqueID());
//  check(float_a.GetUniqueID() + 1 == float_b.GetUniqueID());
//  check(int_a.GetUniqueID() == float_a.GetUniqueID());
//}
