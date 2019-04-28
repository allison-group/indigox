#include "serialise.hpp"
#include "vec_avx.hpp"

#include <cstdint>
#include <cstdlib>
#include <new>

#ifndef INDIGOX_UTILS_ATOMIC_COORDINATES_HPP
#define INDIGOX_UTILS_ATOMIC_COORDINATES_HPP

namespace indigox {

  // Struct for packaging a given set of coordinates
  struct Coordinates {
    double x, y, z;
  };

#define count_per_size 3
#define scale_factor 2
  struct AtomicCoordinates {

    template <typename Archive>
    typename std::enable_if<cereal::traits::is_input_serializable<cereal::BinaryData<double>, Archive>::value
    || cereal::traits::is_output_serializable<cereal::BinaryData<double>, Archive>::value, void>::type
    serialise(Archive& archive, const uint32_t) {
      archive(size);
      if (INDIGOX_IS_INPUT_ARCHIVE(Archive)) Reserve(size);
      archive(cereal::binary_data(data, size * sizeof(double)));
      archive(cereal::binary_data(&data[max_size], size * sizeof(double)));
      archive(cereal::binary_data(&data[2 * max_size], size * sizeof(double)));
    }
    
    template <typename Archive>
    typename std::enable_if<!cereal::traits::is_input_serializable<cereal::BinaryData<double>, Archive>::value
    && !cereal::traits::is_output_serializable<cereal::BinaryData<double>, Archive>::value, void>::type
    serialise(Archive& archive, const uint32_t) {
      archive(INDIGOX_SERIAL_NVP("num", size));
      if (INDIGOX_IS_INPUT_ARCHIVE(Archive)) Reserve(size);
      for (uint32_t i = 0; i < size; ++ i) {
        archive(data[i]);
        archive(data[i + max_size]);
        archive(data[i + 2 * max_size]);
      }
    }
    
    AtomicCoordinates()
        : unaligned_mem(nullptr), data(nullptr), max_size(0), size(0) {}

    AtomicCoordinates(uint32_t count) : AtomicCoordinates() { Reserve(count); }

    AtomicCoordinates(AtomicCoordinates &other) : max_size(0), size(0) {
      Reserve(other.size);
      if (max_size) {
        memcpy(&data[max_size * 0], &other.data[other.max_size * 0], other.size * sizeof(double));
        memcpy(&data[max_size * 1], &other.data[other.max_size * 1], other.size * sizeof(double));
        memcpy(&data[max_size * 2], &other.data[other.max_size * 2], other.size * sizeof(double));
      }
      size = other.size;
    }

    AtomicCoordinates(AtomicCoordinates &&other)
        : unaligned_mem(other.unaligned_mem), data(other.data),
          max_size(other.max_size), size(other.size) {
      other.unaligned_mem = nullptr;
      other.data = nullptr;
    }

    ~AtomicCoordinates() {
      delete[] unaligned_mem;
      max_size = 0;
      size = 0;
    }

    AtomicCoordinates &operator=(AtomicCoordinates &other) {
      // Self check
      if (&other == this) return *this;

      // Release
      delete[] unaligned_mem;
      data = nullptr;

      // Reacquire and copy
      Reserve(other.size);
      if (max_size) {
        memcpy(&data[max_size * 0], &other.data[other.max_size * 0], other.size * sizeof(double));
        memcpy(&data[max_size * 1], &other.data[other.max_size * 1], other.size * sizeof(double));
        memcpy(&data[max_size * 2], &other.data[other.max_size * 2], other.size * sizeof(double));
      }
      size = other.size;

      return *this;
    }

    AtomicCoordinates &operator=(AtomicCoordinates &&other) {
      // Self check
      if (&other == this) return *this;

      // Release
      delete[] unaligned_mem;
      data = nullptr;

      // Steal from other
      unaligned_mem = other.unaligned_mem;
      data = other.data;
      other.unaligned_mem = nullptr;
      other.data = nullptr;

      return *this;
    }

    Coordinates operator[](uint32_t index) {
      return {data[index + 0 * max_size],
              data[index + 1 * max_size],
              data[index + 2 * max_size]};
    }

    double* x_vals() { return data; }
    double* y_vals() { return &data[max_size]; }
    double* z_vals() { return &data[2 * max_size]; }
    
    void SetCoordinates(uint32_t index, double x, double y, double z) {
      data[index + 0 * max_size] = x;
      data[index + 1 * max_size] = y;
      data[index + 2 * max_size] = z;
    }

    // Reserve space for at least count xyz values
    void Reserve(uint32_t count) {
      if (count <= max_size) return;
      AllocateMemory(count);
    }

    // Erase the coordinates at index
    // Replaces the values at index by the values at last position
    // Decrements the size
    void Erase(uint32_t index) {
      // If erasing the last coordinates, just decrement and zero out
      uint32_t last_pos = size - 1;
      if (index != last_pos) {
        data[index + 0 * max_size] = data[last_pos + 0 * max_size];
        data[index + 1 * max_size] = data[last_pos + 1 * max_size];
        data[index + 2 * max_size] = data[last_pos + 2 * max_size];
      }
      data[last_pos + 0 * max_size] = 0.;
      data[last_pos + 1 * max_size] = 0.;
      data[last_pos + 2 * max_size] = 0.;
      --size;
    }

    // Sticks the coordinates at the last position. Expands memory if full.
    // Returns the index of the appended values
    uint32_t Append(double x, double y, double z) {
      if (!data) AllocateMemory(12); // Just start at 12. no real reason.
      if (size == max_size) AllocateMemory(max_size * scale_factor);

      // if append opens up a new group, zero it all out first
      if (!(size % 4)) {
        memset(&data[size], 0, sizeof(double) * 4);
        memset(&data[size + max_size], 0, sizeof(double) * 4);
        memset(&data[size + 2 * max_size], 0, sizeof(double) * 4);
      }

      data[size] = x;
      data[size + max_size] = y;
      data[size + 2 * max_size] = z;

      return size++;
    }

  private:
    char *unaligned_mem = nullptr; // Unaligned raw memory
    char *old_unaligned = nullptr; // Old unaligned raw memory after resizing

  public:
    double *data = nullptr; // Aligned type-d memory

    uint32_t max_size = 0; // Maximum number of xyz can hold
    uint32_t size = 0;     // Current number of held xyz

  private:
    // Allocates space for at least count xyz.
    // Rounds count up to nearest multiple of 4.
    // Copies current memory into new.
    // Puts current into old and leaves for calling method to free
    void AllocateMemory(uint32_t count) {
      const uint32_t align = sizeof(Vec4d) - 1; // align on 32 byte boundary

      // Allocate new memory block
      if (count % 4)
        count += 4 - (count % 4); // Round up to nearest multiple of 4
      uint32_t new_size = count * count_per_size * sizeof(double);
      char *new_mem = new char[new_size + align];
      if (!new_mem) throw std::bad_alloc();
      double *new_data =
          (double *)(((uintptr_t)new_mem + align) & ~(uintptr_t)(align));

      // Copy existing data
      if (data) {
        memcpy(new_data, data, size * sizeof(double));
        memcpy(&new_data[count], &data[max_size], size * sizeof(double));
        memcpy(&new_data[count * 2], &data[max_size * 2], size * sizeof(double));
      }

      // Save old data.
      old_unaligned = unaligned_mem;
      unaligned_mem = new_mem;
      data = new_data;
      max_size = count;

      // Free old data
      if (old_unaligned) {
        delete[] old_unaligned;
        old_unaligned = nullptr;
      }
    }
  };

#undef count_per_size
#undef scale_factor
  
} // namespace indigox

#endif /* INDIGOX_UTILS_ATOMIC_COORDINATES_HPP */
