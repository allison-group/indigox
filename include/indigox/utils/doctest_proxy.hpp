#ifdef INDIGOX_DISABLE_TESTS
#define DOCTEST_CONFIG_DISABLE
#endif

#define DOCTEST_CONFIG_NO_SHORT_MACRO_NAMES
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include <doctest.h>

#define test_case         DOCTEST_TEST_CASE
#define test_case_fixture DOCTEST_TEST_CASE_FIXTURE
#define subcase           DOCTEST_SUBCASE
#define test_suite        DOCTEST_TEST_SUITE
#define test_suite_open   DOCTEST_TEST_SUITE_BEGIN
#define test_suite_close  DOCTEST_TEST_SUITE_END
#define check_throws      DOCTEST_CHECK_THROWS
#define check_throws_as   DOCTEST_CHECK_THROWS_AS
#define check_nothrow     DOCTEST_CHECK_NOTHROW

#define check             DOCTEST_FAST_CHECK_UNARY
#define check_false       DOCTEST_FAST_CHECK_UNARY_FALSE
#define check_eq          DOCTEST_FAST_CHECK_EQ
#define check_ne          DOCTEST_FAST_CHECK_NE
#define check_gt          DOCTEST_FAST_CHECK_GT
#define check_lt          DOCTEST_FAST_CHECK_LT

#ifndef INDIGOX_DISABLE_TESTS
using approximately = doctest::Approx;
using should_fail = doctest::should_fail;
using expected_fails = doctest::expected_failures;
#endif

template <class T1, class T2>
struct TTPair {
  using t1 = T1;
  using t2 = T2;
};

#include "serialise.hpp"
template <class U>
using ixserial = doctest::Types<TTPair<__hr_in,U>, TTPair<__cr_in,U>>;
DOCTEST_TYPE_TO_STRING(__hr_in);
DOCTEST_TYPE_TO_STRING(__cr_in);
