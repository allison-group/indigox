#ifndef INDIGOX_UTILS_NUMERICS_HPP
#define INDIGOX_UTILS_NUMERICS_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

namespace indigox {

  template <typename InputIter, typename T>
  size_t Combinations(InputIter first, InputIter last, size_t r,
                      std::vector<std::vector<T>> &out) {
    out.clear();
    std::vector<T> pool(first, last);
    if (r > pool.size())
      return 0;
    std::vector<size_t> indices(r);
    std::iota(indices.begin(), indices.end(), 0);
    out.emplace_back(std::vector<T>());
    for (size_t i : indices)
      out.back().emplace_back(pool[i]);
    std::vector<size_t> tmp(r);
    std::iota(tmp.begin(), tmp.end(), 0);
    std::reverse(tmp.begin(), tmp.end());
    while (true) {
      size_t i = 0;
      bool broken = false;
      for (size_t x : tmp) {
        if (indices[x] != x + pool.size() - r) {
          i = x;
          broken = true;
          break;
        }
      }
      if (!broken)
        break;
      ++indices[i];
      std::vector<size_t> tmp2(r - i - 1);
      std::iota(tmp2.begin(), tmp2.end(), i + 1);
      for (size_t j : tmp2)
        indices[j] = indices[j - 1] + 1;
      out.emplace_back(std::vector<T>());
      for (size_t i : indices)
        out.back().emplace_back(pool[i]);
    }
    return out.size();
  }

  /*! \brief Calculate the mean of a range of numbers.
   *  \tparam Iter type of the iterator defining the range.
   *  \param begin,end iterators defining the range.
   *  \return the mean of the range. */
  template <class Iter> double CalculateMean(Iter begin, Iter end) {
    double sum = std::accumulate(begin, end, 0.0);
    return sum / std::distance(begin, end);
  }

  /*! \brief Calculate the median of a range of numbers.
   *  \details If the range has an odd length, the median is the middle item in
   *  the sorted range. If it has an even length, the median is the mean of the
   *  two central items. If the range is empty, the median is undefined.
   *  \tparam RandomIter random access itertor defining the range.
   *  \param begin,end iterators defining the range.
   *  \return the median of the range.
   *  \throws std::runtime_error if the range is empty. */
  template <class RandomIter>
  double CalculateMedian(RandomIter begin, RandomIter end) {
    if (begin == end)
      throw std::runtime_error("Median of empty range is undefined");
    size_t sz = end - begin;
    size_t mid = sz / 2;
    RandomIter target = begin + mid;
    std::nth_element(begin, target, end);

    if (sz % 2)
      return static_cast<double>(*target);
    else {
      double x = static_cast<double>(*target);
      RandomIter next = std::max_element(begin, target);
      return (x + static_cast<double>(*next)) / 2.0;
    }
  }

  template <class Iter>
  double CalculateStandardDeviation(Iter begin, Iter end) {
    double mean = CalculateMean(begin, end);
    std::vector<double> d(std::distance(begin, end));
    std::transform(begin, end, d.begin(),
                   [mean](double x) { return x - mean; });
    double sum_sq = std::inner_product(d.begin(), d.end(), d.begin(), 0.0);
    return std::sqrt(sum_sq / d.size());
  }

} // namespace indigox

#endif /* INDIGOX_UTILS_NUMERICS_HPP */
