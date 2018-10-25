#include <utility>

#ifndef INDIGOX_UTILS_TRIPLE_HPP
#define INDIGOX_UTILS_TRIPLE_HPP

#include "serialise.hpp"

namespace stdx {    // extended std namespace
  
  template <class _T1, class _T2, class _T3>
  struct triple
  {
    typedef _T1 first_type;
    typedef _T2 second_type;
    typedef _T3 third_type;
    
    _T1 first;
    _T2 second;
    _T3 third;
    
    triple(triple const&) = default;
    triple(triple&&) = default;
    
    template <bool _Val>
    using _EnableB = typename std::enable_if<_Val, bool>::type;
    
    struct _CheckArgs {
      template <class _U1, class _U2, class _U3>
      static constexpr bool __enable_default() {
        return std::is_default_constructible<_U1>::value
        && std::is_default_constructible<_U2>::value
        && std::is_default_constructible<_U3>::value;
      }
      
      template <class _U1, class _U2, class _U3>
      static constexpr bool __enable_explicit() {
        return std::is_constructible<first_type, _U1>::value
        && std::is_constructible<second_type, _U2>::value
        && std::is_constructible<third_type, _U3>::value
        && (!std::is_convertible<_U1, first_type>::value
            || !std::is_convertible<_U2, second_type>::value
            || !std::is_convertible<_U3, third_type>::value);
      }
      
      template <class _U1, class _U2, class _U3>
      static constexpr bool __enable_implicit() {
        return std::is_constructible<first_type, _U1>::value
        && std::is_constructible<second_type, _U2>::value
        && std::is_constructible<third_type, _U3>::value
        && std::is_convertible<_U1, first_type>::value
        && std::is_convertible<_U2, second_type>::value
        && std::is_convertible<_U3, third_type>::value;
      }
    };
    
    template <bool _MaybeEnable>
    using _CheckArgsDep = typename std::conditional<
    _MaybeEnable, _CheckArgs, std::__check_tuple_constructor_fail>::type;
    
    struct _CheckTupleLikeConstructor {
      template <class _Tuple>
      static constexpr bool __enable_implicit() {
        return std::__tuple_convertible<_Tuple, triple>::value;
      }
      
      template <class _Tuple>
      static constexpr bool __enable_explicit() {
        return std::__tuple_constructible<_Tuple, triple>::value
        && !std::__tuple_convertible<_Tuple, triple>::value;
      }
      
      template <class _Tuple>
      static constexpr bool __enable_assign() {
        return std::__tuple_assignable<_Tuple, triple>::value;
      }
    };
    
    template <class _Tuple>
    using _CheckTLC = typename std::conditional<
    std::__tuple_like_with_size<_Tuple, 3>::value
    && !std::is_same<typename std::decay<_Tuple>::type, triple>::value,
    _CheckTupleLikeConstructor,
    std::__check_tuple_constructor_fail
    >::type;
    
    template<bool _Dummy = true, _EnableB<
    _CheckArgsDep<_Dummy>::template __enable_default<_T1, _T2, _T3>()
    > = false>
    inline constexpr triple() : first(), second(), third() {}
    
    template <bool _Dummy = true, _EnableB<
    _CheckArgsDep<_Dummy>::template
    __enable_explicit<_T1 const&, _T2 const&, _T3 const&>()
    > = false>
    inline constexpr explicit
    triple(_T1 const& __t1, _T2 const& __t2, _T3 const& __t3)
    : first(__t1), second(__t2), third(__t3) {}
    
    template<bool _Dummy = true, _EnableB<
    _CheckArgsDep<_Dummy>::template
    __enable_implicit<_T1 const&, _T2 const&, _T3 const&>()
    > = false>
    inline constexpr
    triple(_T1 const& __t1, _T2 const& __t2, _T3 const& __t3)
    : first(__t1), second(__t2), third(__t3) {}
    
    template<class _U1, class _U2, class _U3, _EnableB<
    _CheckArgs::template __enable_explicit<_U1, _U2, _U3>()
    > = false>
    inline constexpr
    explicit triple(_U1&& __u1, _U2&& __u2, _U3&& __u3)
    : first(std::forward<_U1>(__u1)), second(std::forward<_U2>(__u2)),
    third(std::forward<_U3>(__u3)) {}
    
    template<class _U1, class _U2, class _U3, _EnableB<
    _CheckArgs::template __enable_implicit<_U1, _U2, _U3>()
    > = false>
    inline constexpr
    triple(_U1&& __u1, _U2&& __u2, _U3&& __u3)
    : first(std::forward<_U1>(__u1)), second(std::forward<_U2>(__u2)),
    third(std::forward<_U3>(__u3))  {}
    
    template<class _U1, class _U2, class _U3, _EnableB<
    _CheckArgs::template __enable_explicit<_U1 const&, _U2 const&, _U3 const&>()
    > = false>
    inline constexpr
    explicit triple(triple<_U1, _U2, _U3> const& __p)
    : first(__p.first), second(__p.second), third(__p.third) {}
    
    template<class _U1, class _U2, class _U3, _EnableB<
    _CheckArgs::template __enable_implicit<_U1 const&, _U2 const&, _U3 const&>()
    > = false>
    inline constexpr
    triple(triple<_U1, _U2, _U3> const& __p)
    : first(__p.first), second(__p.second), third(__p.third) {}
    
    template<class _U1, class _U2, class _U3, _EnableB<
    _CheckArgs::template __enable_explicit<_U1, _U2, _U3>()
    > = false>
    inline constexpr
    explicit triple(triple<_U1, _U2, _U3>&&__p)
    : first(std::forward<_U1>(__p.first)), second(std::forward<_U2>(__p.second)),
    third(std::forward<_U3>(__p.third)) {}
    
    template<class _U1, class _U2, class _U3, _EnableB<
    _CheckArgs::template __enable_implicit<_U1, _U2, _U3>()
    > = false>
    inline constexpr
    triple(triple<_U1, _U2, _U3>&& __p)
    : first(std::forward<_U1>(__p.first)), second(std::forward<_U2>(__p.second)),
    third(std::forward<_U3>(__p.third)) {}
    
    inline
    triple& operator=(typename std::conditional<
                      std::is_copy_assignable<first_type>::value &&
                      std::is_copy_assignable<second_type>::value &&
                      std::is_copy_assignable<third_type>::value,
                      triple, std::__nat>::type const& __p)
    noexcept(std::is_nothrow_copy_assignable<first_type>::value &&
             std::is_nothrow_copy_assignable<second_type>::value &&
             std::is_nothrow_copy_assignable<third_type>::value)
    {
      first = __p.first;
      second = __p.second;
      third = __p.third;
      return *this;
    }
    
    inline
    triple& operator=(typename std::conditional<
                      std::is_move_assignable<first_type>::value &&
                      std::is_move_assignable<second_type>::value &&
                      std::is_move_assignable<third_type>::value,
                      triple, std::__nat>::type&& __p)
    noexcept(std::is_nothrow_move_assignable<first_type>::value &&
             std::is_nothrow_move_assignable<second_type>::value &&
             std::is_nothrow_move_assignable<third_type>::value)
    {
      first = std::forward<first_type>(__p.first);
      second = std::forward<second_type>(__p.second);
      third = std::forward<third_type>(__p.third);
      return *this;
    }
    
    inline
    void
    swap(triple& __p) noexcept(std::__is_nothrow_swappable<first_type>::value &&
                               std::__is_nothrow_swappable<second_type>::value &&
                               std::__is_nothrow_swappable<third_type>::value)
    {
      using std::swap;
      swap(first,  __p.first);
      swap(second, __p.second);
      swap(third, __p.third);
    }
    
    template <typename Archive>
    void serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("first", first),
              INDIGOX_SERIAL_NVP("second", second),
              INDIGOX_SERIAL_NVP("third", third));
    }
  };
  
  template <class _T1, class _T2, class _T3>
  inline constexpr
  bool
  operator==(const triple<_T1,_T2,_T3>& __x, const triple<_T1,_T2,_T3>& __y)
  {
    return __x.first == __y.first &&
           __x.second == __y.second &&
           __x.third == __y.third;
  }
  
  template <class _T1, class _T2, class _T3>
  inline constexpr
  bool
  operator!=(const triple<_T1,_T2,_T3>& __x, const triple<_T1,_T2,_T3>& __y)
  {
    return !(__x == __y);
  }
  
  template <class _T1, class _T2, class _T3>
  inline constexpr
  bool
  operator< (const triple<_T1,_T2,_T3>& __x, const triple<_T1,_T2,_T3>& __y)
  {
    return __x.first < __y.first ||
    (!(__y.first < __x.first) && __x.second < __y.second) ||
    (!(__y.first < __x.first) && !(__y.second < __x.second) && __x.third < __y.third);
  }
  
  template <class _T1, class _T2, class _T3>
  inline constexpr
  bool
  operator> (const triple<_T1,_T2,_T3>& __x, const triple<_T1,_T2,_T3>& __y)
  {
    return __y < __x;
  }
  
  template <class _T1, class _T2, class _T3>
  inline constexpr
  bool
  operator>=(const triple<_T1,_T2,_T3>& __x, const triple<_T1,_T2,_T3>& __y)
  {
    return !(__x < __y);
  }
  
  template <class _T1, class _T2, class _T3>
  inline constexpr
  bool
  operator<=(const triple<_T1,_T2,_T3>& __x, const triple<_T1,_T2,_T3>& __y)
  {
    return !(__y < __x);
  }
  
  template <class _T1, class _T2, class _T3>
  inline
  typename std::enable_if
  <
  std::__is_swappable<_T1>::value &&
  std::__is_swappable<_T2>::value &&
  std::__is_swappable<_T3>::value,
  void
  >::type
  swap(triple<_T1, _T2, _T3>& __x, triple<_T1, _T2, _T3>& __y)
  noexcept((std::__is_nothrow_swappable<_T1>::value &&
            std::__is_nothrow_swappable<_T2>::value &&
            std::__is_nothrow_swappable<_T3>::value))
  {
    __x.swap(__y);
  }
  
  template <class _Tp>
  struct __make_triple_return_impl
  {
    typedef _Tp type;
  };
  
  template <class _Tp>
  struct __make_triple_return_impl<std::reference_wrapper<_Tp>>
  {
    typedef _Tp& type;
  };
  
  template <class _Tp>
  struct __make_triple_return
  {
    typedef typename __make_triple_return_impl<typename std::decay<_Tp>::type>::type type;
  };
  
  template <class _T1, class _T2, class _T3>
  inline constexpr
  triple<typename __make_triple_return<_T1>::type,
         typename __make_triple_return<_T2>::type,
         typename __make_triple_return<_T3>::type>
  make_triple(_T1&& __t1, _T2&& __t2, _T3&& __t3)
  {
    return triple<typename __make_triple_return<_T1>::type,
                  typename __make_triple_return<_T2>::type,
                  typename __make_triple_return<_T3>::type>
    (std::forward<_T1>(__t1), std::forward<_T2>(__t2), std::forward<_T3>(__t3));
  }
}

#endif /* INDIGOX_UTILS_TRIPLE_HPP */
