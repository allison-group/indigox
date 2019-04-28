#include <utility>
#include <x86intrin.h>

#ifndef INDIGOX_UTILS_VEC_AVX_HPP
#define INDIGOX_UTILS_VEC_AVX_HPP

namespace indigox {

  // AVX vector of 4 double floating point values
  struct Vec4d {
    using real_t = double;
    using vec_t = __m256d;

    vec_t raw;

    static Vec4d ZeroVector() { return _mm256_setzero_pd(); }

    // Default constructor
    Vec4d() {}

    // Construct with given value in all slots
    Vec4d(real_t val) : raw(_mm256_set1_pd(val)) {}

    // Construct with given values
    Vec4d(real_t v0, real_t v1, real_t v2, real_t v3)
        : raw(_mm256_setr_pd(v0, v1, v2, v3)) {}

    // Construct from pointer. Assume unalligned
    Vec4d(real_t *arr) : raw(_mm256_loadu_pd(arr)) {}

    // Convert from raw type
    Vec4d(vec_t r) : raw(r) {}

    // implicit convert to raw type
    operator vec_t() { return raw; }

    // Copy construct
    Vec4d(const Vec4d &other) : raw(other.raw) {}

    // Move construct
    Vec4d(Vec4d &&other) : raw(std::move(other.raw)) {}

    // Copy assign
    Vec4d &operator=(const Vec4d &other) {
      if (&other != this) raw = other.raw;
      return *this;
    }

    // Move assign
    Vec4d &operator=(Vec4d &&other) {
      if (&other != this) raw = std::move(other.raw);
      return *this;
    }

    // Assign from vec_T
    Vec4d &operator=(const vec_t &other) {
      raw = other;
      return *this;
    }

    // Assign from single real_t
    Vec4d &operator=(real_t other) {
      raw = _mm256_set1_pd(other);
      return *this;
    }

    // Assign from array of real_t. Assume unaligned
    Vec4d &operator=(real_t *arr) {
      raw = _mm256_loadu_pd(arr);
      return *this;
    }

    // Store value into an array
    void save(real_t *arr) { _mm256_storeu_pd(arr, raw); }

    // Get the i-th value from the vector. Performs a range breakdown
    real_t operator[](uint32_t i) { return raw[i & 7]; }

    // Add a Vec4d to this and return a new Vec4d
    inline Vec4d operator+(const Vec4d &other) {
      return _mm256_add_pd(raw, other.raw);
    }

    // Add a scalar to this and return new Vec4d
    inline Vec4d operator+(real_t other) {
      return _mm256_add_pd(raw, Vec4d(other));
    }

    // Add in place a Vec4d to this
    inline Vec4d &operator+=(const Vec4d &other) {
      raw = _mm256_add_pd(raw, other.raw);
      return *this;
    }

    // Add in place a scalar to this
    inline Vec4d &operator+=(real_t other) {
      raw = _mm256_add_pd(raw, Vec4d(other));
      return *this;
    }

    // Subtract a Vec4d from this and return new Vec4d
    inline Vec4d operator-(const Vec4d &other) {
      return _mm256_sub_pd(raw, other.raw);
    }

    // Subtract a scalar from this and return new Vec4d
    inline Vec4d operator-(real_t other) {
      return _mm256_sub_pd(raw, Vec4d(other));
    }

    // Subtract in place a Vec4d from this
    inline Vec4d &operator-=(const Vec4d &other) {
      raw = _mm256_sub_pd(raw, other.raw);
      return *this;
    }

    // Subtract in place a scalar from this
    inline Vec4d &operator-=(real_t other) {
      raw = _mm256_sub_pd(raw, Vec4d(other));
      return *this;
    }

    // Multiply by a Vec4d (element by element)
    inline Vec4d operator*(const Vec4d &other) {
      return _mm256_mul_pd(raw, other.raw);
    }

    // Multiply each element by a scalar
    inline Vec4d operator*(real_t other) {
      return _mm256_mul_pd(raw, Vec4d(other));
    }

    // Multiply in place a Vec4d element by element
    inline Vec4d &operator*=(const Vec4d &other) {
      raw = _mm256_mul_pd(raw, other.raw);
      return *this;
    }

    // Multiply in place each element by a scalar
    inline Vec4d &operator*=(real_t other) {
      raw = _mm256_mul_pd(raw, Vec4d(other));
      return *this;
    }
    
    // Divide by a Vec4d (element by element)
    inline Vec4d operator/(const Vec4d &other) {
      return _mm256_div_pd(raw, other.raw);
    }
    
    //
  };

  /****************************
   *    Operators for Vec4d    *
   ****************************/

  // Out of class add a Vec4d to a scalar return new Vec4d
  inline Vec4d operator+(Vec4d::real_t a, const Vec4d &b) {
    return _mm256_add_pd(Vec4d(a), b.raw);
  }

  // Out of class subtract a Vec4d from a scalr return new Vec4d
  inline Vec4d operator-(Vec4d::real_t a, const Vec4d &b) {
    return _mm256_sub_pd(Vec4d(a), b.raw);
  }

  // Out of class multiply a Vec4d by a scalar
  inline Vec4d operator*(Vec4d::real_t a, const Vec4d &b) {
    return _mm256_mul_pd(Vec4d(a), b.raw);
  }

} // namespace indigox

#endif /* INDIGOX_UTILS_VEC3_HPP */
