#include <type_traits>
#include <utility>

#ifndef INDIGOX_UTILS_QUAD_HPP
#define INDIGOX_UTILS_QUAD_HPP

#include "serialise.hpp"

namespace stdx {    // extended std namespace
  
  template <class _T1, class _T2=_T1, class _T3=_T2, class _T4=_T3>
  struct quad
  {
    typedef _T1 first_type;
    typedef _T2 second_type;
    typedef _T3 third_type;
    typedef _T4 fourth_type;
    
    _T1 first;
    _T2 second;
    _T3 third;
    _T4 fourth;
    
    quad(quad const&) = default;
    quad(quad&&) = default;
    
    inline constexpr quad() : first(), second(), third(), fourth() {}
    
    inline constexpr
    quad(_T1 const& __t1, _T2 const& __t2, _T3 const& __t3, _T4 const& __t4)
    : first(__t1), second(__t2), third(__t3) , fourth(__t4)  {}
    
    template<class _U1, class _U2, class _U3, class _U4>
    inline constexpr
    quad(_U1&& __u1, _U2&& __u2, _U3&& __u3, _U4&& __u4)
    : first(std::forward<_U1>(__u1)), second(std::forward<_U2>(__u2)),
    third(std::forward<_U3>(__u3)), fourth(std::forward<_U4>(__u4))  {}
    
    template<class _U1, class _U2, class _U3, class _U4>
    inline constexpr
    quad(quad<_U1, _U2, _U3, _U4> const& __p)
    : first(__p.first), second(__p.second), third(__p.third), fourth(__p.fourth) {}
    
    template<class _U1, class _U2, class _U3, class _U4>
    inline constexpr
    quad(quad<_U1, _U2, _U3, _U4>&& __p)
    : first(std::forward<_U1>(__p.first)), second(std::forward<_U2>(__p.second)),
    third(std::forward<_U3>(__p.third)), fourth(std::forward<_U4>(__p.fourth)) {}
    
    inline
    quad& operator=(typename std::conditional<
                    std::is_copy_assignable<first_type>::value &&
                    std::is_copy_assignable<second_type>::value &&
                    std::is_copy_assignable<third_type>::value &&
                    std::is_copy_assignable<fourth_type>::value,
                    quad, stdx::typeless>::type const& __p)
    noexcept(std::is_nothrow_copy_assignable<first_type>::value &&
             std::is_nothrow_copy_assignable<second_type>::value &&
             std::is_nothrow_copy_assignable<third_type>::value &&
             std::is_nothrow_copy_assignable<fourth_type>::value)
    {
      first = __p.first;
      second = __p.second;
      third = __p.third;
      fourth = __p.fourth;
      return *this;
    }
    
    inline
    quad& operator=(typename std::conditional<
                    std::is_move_assignable<first_type>::value &&
                    std::is_move_assignable<second_type>::value &&
                    std::is_move_assignable<third_type>::value &&
                    std::is_move_assignable<fourth_type>::value,
                    quad, stdx::typeless>::type&& __p)
    noexcept(std::is_nothrow_move_assignable<first_type>::value &&
             std::is_nothrow_move_assignable<second_type>::value &&
             std::is_nothrow_move_assignable<third_type>::value &&
             std::is_nothrow_move_assignable<fourth_type>::value)
    {
      first = std::forward<first_type>(__p.first);
      second = std::forward<second_type>(__p.second);
      third = std::forward<third_type>(__p.third);
      fourth = std::forward<fourth_type>(__p.fourth);
      return *this;
    }
    
    inline
    void
    swap(quad& __p) noexcept(std::is_nothrow_swappable<first_type>::value &&
                             std::is_nothrow_swappable<second_type>::value &&
                             std::is_nothrow_swappable<third_type>::value&&
                             std::is_nothrow_swappable<fourth_type>::value)
    {
      using std::swap;
      swap(first,  __p.first);
      swap(second, __p.second);
      swap(third, __p.third);
      swap(fourth, __p.fourth);
    }
    
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("first", first),
              INDIGOX_SERIAL_NVP("second", second),
              INDIGOX_SERIAL_NVP("third", third),
              INDIGOX_SERIAL_NVP("fourth", fourth));
    }
    
  };
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline constexpr
  bool
  operator==(const quad<_T1,_T2,_T3,_T4>& __x, const quad<_T1,_T2,_T3,_T4>& __y)
  {
    return __x.first == __y.first &&
           __x.second == __y.second &&
           __x.third == __y.third &&
           __x.fourth == __y.fourth;
  }
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline constexpr
  bool
  operator!=(const quad<_T1,_T2,_T3,_T4>& __x, const quad<_T1,_T2,_T3,_T4>& __y)
  {
    return !(__x == __y);
  }
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline constexpr
  bool
  operator< (const quad<_T1,_T2,_T3,_T4>& __x, const quad<_T1,_T2,_T3,_T4>& __y)
  {
    return __x.first < __y.first ||
    (!(__y.first < __x.first) && __x.second < __y.second) ||
    (!(__y.first < __x.first) && !(__y.second < __x.second) && __x.third < __y.third) ||
    (!(__y.first < __x.first) && !(__y.second < __x.second) && !(__y.third < __x.third) && __x.fourth < __y.fourth);
  }
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline constexpr
  bool
  operator> (const quad<_T1,_T2,_T3,_T4>& __x, const quad<_T1,_T2,_T3,_T4>& __y)
  {
    return __y < __x;
  }
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline constexpr
  bool
  operator>=(const quad<_T1,_T2,_T3,_T4>& __x, const quad<_T1,_T2,_T3,_T4>& __y)
  {
    return !(__x < __y);
  }
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline constexpr
  bool
  operator<=(const quad<_T1,_T2,_T3,_T4>& __x, const quad<_T1,_T2,_T3,_T4>& __y)
  {
    return !(__y < __x);
  }
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline
  typename std::enable_if
  <
  std::is_swappable<_T1>::value &&
  std::is_swappable<_T2>::value &&
  std::is_swappable<_T3>::value &&
  std::is_swappable<_T4>::value,
  void
  >::type
  swap(quad<_T1, _T2, _T3, _T4>& __x, quad<_T1, _T2, _T3, _T4>& __y)
  noexcept((std::is_nothrow_swappable<_T1>::value &&
            std::is_nothrow_swappable<_T2>::value &&
            std::is_nothrow_swappable<_T3>::value &&
            std::is_nothrow_swappable<_T4>::value))
  {
    __x.swap(__y);
  }
  
  template <class _Tp>
  struct __make_quad_return_impl
  {
    typedef _Tp type;
  };
  
  template <class _Tp>
  struct __make_quad_return_impl<std::reference_wrapper<_Tp>>
  {
    typedef _Tp& type;
  };
  
  template <class _Tp>
  struct __make_quad_return
  {
    typedef typename __make_quad_return_impl<typename std::decay<_Tp>::type>::type type;
  };
  
  template <class _T1, class _T2, class _T3, class _T4>
  inline constexpr
  quad<typename __make_quad_return<_T1>::type,
  typename __make_quad_return<_T2>::type,
  typename __make_quad_return<_T3>::type,
  typename __make_quad_return<_T4>::type>
  make_quad(_T1&& __t1, _T2&& __t2, _T3&& __t3, _T4&& __t4)
  {
    return quad<typename __make_quad_return<_T1>::type,
    typename __make_quad_return<_T2>::type,
    typename __make_quad_return<_T3>::type,
    typename __make_quad_return<_T4>::type>
    (std::forward<_T1>(__t1), std::forward<_T2>(__t2),
     std::forward<_T3>(__t3), std::forward<_T4>(__t4));
  }
}

#endif /* INDIGOX_UTILS_TRIPLE_HPP */
