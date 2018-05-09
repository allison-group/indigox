#include <indigox/classes/periodictable.hpp>
#include <indigox/utils/helpers.hpp>

namespace indigox {
  
  PeriodicTable GetPeriodicTable() {
    static PeriodicTable instance = PeriodicTable();
    if (!instance) {
      instance.reset(new IXPeriodicTable());
      instance->GeneratePeriodicTable();
    }
    return instance;
  }
}
