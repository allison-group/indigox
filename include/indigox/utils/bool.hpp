#ifndef INDIGOX_UTILS_BOOL_HPP
#define INDIGOX_UTILS_BOOL_HPP

#include <string>
namespace indigox::utils {
  struct IXBool {
    const bool success;
    const std::string reason;
    IXBool(bool p, std::string r) : success(p), reason(r) { }
    operator bool() const { return success; }
  };
}

#endif /* INDIGOX_UTILS_BOOL_HPP */
