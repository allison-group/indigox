/** @file filereader.hpp
 *  @brief Declaration of the FileReader class.
 *  @author Ivan Welsh
 *  @date 13 December 2017
 *  @lastmodify 6 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#ifndef INDIGOX_UTILS_FILEREADER_HPP
#define INDIGOX_UTILS_FILEREADER_HPP

#include <vector>

#include "numerics.hpp"

namespace indigox {
  namespace utils {
    
    /** @class FileReader filereader.hpp
     *  @brief Class for reading simple text file from disk.
     *  @details Loads a simple text file from disk.
     *  @since 0.1
     */
    class FileReader {
      
    public:
      FileReader() = delete;  // No default constructor
      /// @brief Normal constructor.
      FileReader(const string_&);
      
      /// @brief Reads the file from disk.
      void GetAllLines(std::vector<string_>&);
      inline void SetFilePath(const string_& p) { path_ = p; }
      
    private:
      string_ path_;
      
    };
    
  }
}

#endif /* INDIGOX_UTILS_FILEREADER_HPP */
