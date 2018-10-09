#ifndef INDIGOX_UTILS_NUMERICS_HPP
#define INDIGOX_UTILS_NUMERICS_HPP

#include <cstdint>
#include <string>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <cmath>

namespace indigox {
//  using char_ = int8_t;
//  using uchar_ = uint8_t;
//  using short_ = int16_t;
//  using ushort_ = uint16_t;
//  using int_ = int32_t;
//  using uint_ = uint32_t;
//  using long_ = int64_t;
//  using ulong_ = uint64_t;
//  using size_ = size_t;
//  using uid_ = uint64_t;
  
  /*! \brief Calculate the mean of a range of numbers.
   *  \tparam Iter type of the iterator defining the range.
   *  \param begin,end iterators defining the range.
   *  \return the mean of the range. */
  template <class Iter>
  double CalculateMean(Iter begin, Iter end) {
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
    if (begin == end) throw std::runtime_error("Median of empty range is undefined");
    size_t sz = end - begin;
    size_t mid = sz / 2;
    RandomIter target = begin + mid;
    std::nth_element(begin, target, end);
    
    if (sz % 2) return static_cast<double>(*target);
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
    std::transform(begin, end, d.begin(), [mean](double x){return x - mean;});
    double sum_sq = std::inner_product(d.begin(), d.end(), d.begin(), 0.0);
    return std::sqrt(sum_sq / d.size());
  }
  
  
} // namespace indigox

#endif /* INDIGOX_UTILS_NUMERICS_HPP */
