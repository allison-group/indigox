#ifndef INDIGOX_UTILS_HELPERS_HPP
#define INDIGOX_UTILS_HELPERS_HPP

namespace indigox {
  class IXPeriodicTable;
  typedef std::shared_ptr<IXPeriodicTable> PeriodicTable;
  
  PeriodicTable GetPeriodicTable();
}

#endif /* INDIGOX_UTILS_HELPERS_HPP */
