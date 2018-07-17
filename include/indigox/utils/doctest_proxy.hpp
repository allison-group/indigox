#ifndef doctest_proxy_h
#define doctest_proxy_h

#define DOCTEST_CONFIG_NO_SHORT_MACRO_NAMES
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include <doctest.h>

#define test_case       DOCTEST_TEST_CASE
#define subcase         DOCTEST_SUBCASE
#define test_suite      DOCTEST_TEST_SUITE
#define check_throws    DOCTEST_CHECK_THROWS
#define check_throws_as DOCTEST_CHECK_THROWS_AS
#define check_nothrow   DOCTEST_CHECK_NOTHROW

#define check_eq        DOCTEST_FAST_CHECK_EQ
#define check_ne        DOCTEST_FAST_CHECK_NE
#define check_gt        DOCTEST_FAST_CHECK_GT
#define check_lt        DOCTEST_FAST_CHECK_LT
#define check           DOCTEST_FAST_CHECK_UNARY
#define check_not       DOCTEST_FAST_CHECK_UNARY_FALSE

using approximately = doctest::Approx;

#endif /* doctest_proxy_h */
