/** @file filereader.cpp
 *  @brief Implementation of the FileReader class.
 *  @author Ivan Welsh
 *  @date 13 December 2017
 *  @lastmodify 6 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <indigox/utils/filereader.hpp>


namespace indigox::utils {
  
  /** @details Sets the path of a file to read.
   *  @param file path to the file to read.
   */
  FileReader::FileReader(const std::string& file)
  : path_(file) {}
  
  /** @details Loads the given file and reads it itemwise (white space
   *  seperator). Lines begining with the comment character (\#) are
   *  ignored. If the comment character appears in a line, the remainder
   *  of the line is ignored.
   *  @param[out] out_items a vector to store the items from the file.
   */
  void FileReader::GetAllLines(std::vector<std::string> &out_items)
  {
    out_items.clear();
    std::string next_item, line;
    std::ifstream infile;
    infile.open(path_);
    if (!infile) throw std::invalid_argument("File can't be opened.");
    while (std::getline(infile, line)) {
      boost::trim(line);
      out_items.push_back(line);
    }
    infile.close();
  }
  
}
