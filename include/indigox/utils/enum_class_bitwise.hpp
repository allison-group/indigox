#include <type_traits>

#ifndef INDIGOX_UTILS_ENUM_CLASS_BITWISE_HPP
#define INDIGOX_UTILS_ENUM_CLASS_BITWISE_HPP

#define BITWISE_OPERATORS(enum_type)                                           \
  inline enum_type operator|(enum_type l, enum_type r) {                       \
    using T = std::underlying_type_t<enum_type>;                               \
    return static_cast<enum_type>(static_cast<T>(l) | static_cast<T>(r));      \
  }                                                                            \
  inline enum_type &operator|=(enum_type &l, enum_type r) {                    \
    using T = std::underlying_type_t<enum_type>;                               \
    l = static_cast<enum_type>(static_cast<T>(l) | static_cast<T>(r));         \
    return l;                                                                  \
  }                                                                            \
  inline enum_type operator&(enum_type l, enum_type r) {                       \
    using T = std::underlying_type_t<enum_type>;                               \
    return static_cast<enum_type>(static_cast<T>(l) & static_cast<T>(r));      \
  }                                                                            \
  inline enum_type &operator&=(enum_type &l, enum_type r) {                    \
    using T = std::underlying_type_t<enum_type>;                               \
    l = static_cast<enum_type>(static_cast<T>(l) & static_cast<T>(r));         \
    return l;                                                                  \
  }

#endif /* INDIGOX_UTILS_ENUM_CLASS_BITWISE_HPP */
