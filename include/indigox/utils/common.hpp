/** @file common.hpp
 *  @brief Utility functions declarations.
 *  @author Ivan Welsh
 *  @date 21 August 2017
 *  @lastmodify 6 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#ifndef INDIGOX_UTILS_COMMON_HPP
#define INDIGOX_UTILS_COMMON_HPP

#include <string>

namespace indigox::utils {
  
  /// @brief Convert a string to upper case.
  std::string toUpper(const std::string*);
  
  /// @brief Convert a string to lower case.
  std::string toLower(const std::string*);
  
  /// @brief Convert a string to lower case with leading upper case.
  std::string toUpperFirst(const std::string*);
  
  /// @brief Generate a random string.
  std::string randomString(size_t,
                           const std::string chrs = "qwertyuiopasdfghjklzxcvbnmZAQXSWCDEVFRBGTNHYMJUKILOP");
}  // namespace indigox


#endif /* INDIGOX_UTILS_COMMON_HPP */
