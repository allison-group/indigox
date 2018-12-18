#include "triple.hpp"

#include <algorithm>
#include <iterator>
#include <vector>

#ifndef INDIGOX_UTIL_COMBINATRONICS_HPP
#define INDIGOX_UTIL_COMBINATRONICS_HPP

template <class T, class ICT = std::vector<T>> struct CartesianProduct {
  using type = T;
  using innerType = ICT;

  using innerIter = typename innerType::const_iterator;
  using innerIters = stdx::triple<innerIter, innerIter, innerIter>;
  using innerItersC = std::vector<innerIters>;

  innerItersC iters;
  bool finished;

  CartesianProduct() = delete;

  template <class outerIter>
  CartesianProduct(outerIter begin, outerIter end) : finished(false) {
    for (; begin != end; ++begin)
      iters.emplace_back(begin->begin(), begin->end(), begin->begin());
    if (iters.empty())
      finished = true;
  }

  CartesianProduct(innerType &a, innerType &b) : finished(false) {
    iters.emplace_back(a.begin(), a.end(), a.begin());
    iters.emplace_back(b.begin(), b.end(), b.begin());
    if (a.begin() == a.end() || b.begin() == b.end())
      finished = true;
  }

  // Returns true if c contains a product, false otherwise
  bool operator()(innerType &c) {
    c.clear();
    if (finished)
      return false;
    c.reserve(iters.size());
    for (auto it : iters)
      c.push_back(*(it.third));

    for (auto it = iters.begin();;) {
      ++(it->third);
      if (it->third == it->second) {
        if (it + 1 == iters.end())
          finished = true;
        else {
          it->third = it->first;
          ++it;
        }
      } else
        break;
    }
    return true;
  }
};

template <class RandomIt> struct RegionalPermutation {
  // first, second defines range, third defines current
  using IterPair = stdx::triple<RandomIt, RandomIt, RandomIt>;
  using Iters = std::vector<IterPair>;

  RandomIt b, e;
  bool finished, first_pass;
  Iters regions;

  RegionalPermutation() = delete;
  RegionalPermutation(RandomIt begin, RandomIt end)
      : b(begin), e(end), finished(false), first_pass(false) {
  }

  bool AddRegion(RandomIt first, RandomIt last) {
    if (first_pass)
      return false;
    if (first == last)
      return false;
    RandomIt tmp = last;
    if (first == --tmp)
      return false;
    // Check region doesn't overlap an existing region
    for (IterPair &be : regions) {
      if (first < be.first && last <= be.first)
        continue;
      if (first >= be.second && last > be.second)
        continue;
      return false;
    }
    regions.emplace_back(first, last, first);
    return true;
  }

  // Applies the next permutation
  bool operator()() {
    if (!first_pass) {
      first_pass = true;
      for (IterPair &be : regions)
        std::sort(be.first, be.second);
      return true;
    }
    for (IterPair &be : regions) {
      bool fin_region = !std::next_permutation(be.first, be.second);
      if (!fin_region)
        return true;
    }
    return false;
  }
};

#endif /* INDIGOX_UTIL_COMBINATRONICS_HPP */
