#include <EASTL/vector_map.h>

#ifndef INDIGOX_UTILS_SIMPLE_BIMAP_HPP
#define INDIGOX_UTILS_SIMPLE_BIMAP_HPP

namespace indigox::utils {
  template <class L, class R>
  struct SimpleBiMap {
    eastl::vector_map<L, R> left;
    eastl::vector_map<R, L> right;
    
    inline void insert(L l, R r) {
      left.emplace(l, r);
      right.emplace(r, l);
    }
    inline void erase(L l) {
      right.erase(left.at(l));
      left.erase(l);
    }
    inline void erase(R r) {
      left.erase(right.at(r));
      right.erase(r);
    }
    inline void clear() {
      left.clear();
      right.clear();
    }
  };
}

#endif /* INDIGOX_UTILS_SIMPLE_BIMAP_HPP */
