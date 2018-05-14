#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>

#include <indigox/utils/common.hpp>
#include <indigox/utils/numerics.hpp>

using namespace indigox;
using namespace indigox::utils;

struct DummyClassWeakContains {};

BOOST_AUTO_TEST_SUITE(common);

BOOST_AUTO_TEST_CASE(to_upper) {
  string_ test1 = "teststring";
  string_ test2 = "TESTSTRING";
  string_ test3 = "tEsTsTrInG";
  string_ test4 = "1234567890_-+=(){}[]";
  string_ test5 = "1234567890_-+=(){}[]ABCdef";
  BOOST_TEST(ToUpper(test1) == "TESTSTRING");
  BOOST_TEST(ToUpper(test2) == "TESTSTRING");
  BOOST_TEST(ToUpper(test3) == "TESTSTRING");
  BOOST_TEST(ToUpper(test4) == "1234567890_-+=(){}[]");
  BOOST_TEST(ToUpper(test5) == "1234567890_-+=(){}[]ABCDEF");
}

BOOST_AUTO_TEST_CASE(to_lower) {
  string_ test1 = "teststring";
  string_ test2 = "TESTSTRING";
  string_ test3 = "tEsTsTrInG";
  string_ test4 = "1234567890_-+=(){}[]";
  string_ test5 = "1234567890_-+=(){}[]ABCdef";
  BOOST_TEST(ToLower(test1) == "teststring");
  BOOST_TEST(ToLower(test2) == "teststring");
  BOOST_TEST(ToLower(test3) == "teststring");
  BOOST_TEST(ToLower(test4) == "1234567890_-+=(){}[]");
  BOOST_TEST(ToLower(test5) == "1234567890_-+=(){}[]abcdef");
}

BOOST_AUTO_TEST_CASE(to_upper_first) {
  string_ test1 = "teststring";
  string_ test2 = "TESTSTRING";
  string_ test3 = "tEsTsTrInG";
  string_ test4 = "1234567890_-+=(){}[]";
  string_ test5 = "1234567890_-+=(){}[]ABCdef";
  string_ test6 = "test String";
  BOOST_TEST(ToUpperFirst(test1) == "Teststring");
  BOOST_TEST(ToUpperFirst(test2) == "Teststring");
  BOOST_TEST(ToUpperFirst(test3) == "Teststring");
  BOOST_TEST(ToUpperFirst(test4) == "1234567890_-+=(){}[]");
  BOOST_TEST(ToUpperFirst(test5) == "1234567890_-+=(){}[]abcdef");
  BOOST_TEST(ToUpperFirst(test6) == "Test string");
}

BOOST_AUTO_TEST_CASE(random_string) {
  BOOST_TEST(GetRandomString(5, 123456789) == "QOmZX");
  BOOST_TEST(GetRandomString(12) == "MiFhzJxCVKgl");
  BOOST_TEST(GetRandomString(1) == "F");
  BOOST_TEST(GetRandomString(7) == "WMcYLPH");
}

BOOST_AUTO_TEST_CASE(weak_contains) {
  typedef std::shared_ptr<DummyClassWeakContains> shared;
  typedef std::weak_ptr<DummyClassWeakContains> weak;
  std::vector<weak> cases;
  cases.reserve(5);
  shared c1 = shared(new DummyClassWeakContains());
  shared c2 = shared(new DummyClassWeakContains());
  shared c3 = shared(new DummyClassWeakContains());
  shared c4 = shared(new DummyClassWeakContains());
  shared c5 = shared(new DummyClassWeakContains());
  cases.emplace_back(c1); cases.emplace_back(c2); cases.emplace_back(c3);
  cases.emplace_back(c4);
  // Check for existing ones
  BOOST_CHECK(WeakContainsShared(cases.begin(), cases.end(), c3) != cases.end());
  BOOST_CHECK(WeakContainsShared(cases.begin(), cases.end(), c1) == cases.begin());
  // Check for non-existent ones
  BOOST_CHECK(WeakContainsShared(cases.begin(), cases.end(), c5) == cases.end());
  // Check for empty
  BOOST_CHECK(WeakContainsShared(cases.begin(), cases.end(), shared()) == cases.end());
}

BOOST_AUTO_TEST_SUITE_END();
