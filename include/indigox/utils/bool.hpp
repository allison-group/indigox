#ifndef INDIGOX_UTILS_BOOL_HPP
#define INDIGOX_UTILS_BOOL_HPP

#include <string>

namespace indigox::utils {
  /*! \brief A boolean class for indicating a reason for the boolean state. */
  struct IXBool {
    //! \brief The boolean state.
    const bool success;
    //! \brief The reason for the state.
    const std::string reason;
    //! \brief Default constructor.
    IXBool(bool p, std::string r) : success(p), reason(r) { }
    //! \brief Implicit conversion to the primitive bool type.
    operator bool() const { return success; }
  };
}

#endif /* INDIGOX_UTILS_BOOL_HPP */
