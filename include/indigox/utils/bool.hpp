#ifndef INDIGOX_UTILS_BOOL_HPP
#define INDIGOX_UTILS_BOOL_HPP

#include <string>

namespace indigox::utils {
  /*! \brief A boolean class for indicating a reason for the boolean state. */
  struct IXBool {
    enum class Reason {
      // General fails
      BAD_REQUEST,
      // Graph fails
      NO_DUPLICATES,
      NO_SELF_LOOPS,
      NO_PARALLEL_EDGES,
      MISSING_EDGE,
      // Pass
      SUCCESSFUL
    };
    
    //! \brief The boolean state.
    bool success;
    //! \brief The reason for the state.
    Reason reason;
    //! \brief Successful constructor
    IXBool() : success(true), reason(Reason::SUCCESSFUL) { }
    //! \brief Failed constructor.
    IXBool(Reason r) : success(false), reason(r) { }
    //! \brief Implicit conversion to the primitive bool type.
    operator bool() const { return success; }
    //! \brief Assignment operator.
    IXBool& operator=(const IXBool& o) {
      if (&o == this) return *this;
      success = o.success; reason = o.reason;
      return *this;
    }
    //! \brief Compare an IXBool with a Reason
    inline friend bool operator==(const IXBool& b, IXBool::Reason r) {
      return b.reason == r;
    }
    //! \brief Compare a Reason with an IXBool
    inline friend bool operator==(IXBool::Reason r, const IXBool& b) {
      return b.reason == r;
    }
  };
}

#endif /* INDIGOX_UTILS_BOOL_HPP */
