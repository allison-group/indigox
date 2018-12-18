#include <indigox/utils/common.hpp>

#include <algorithm>
#include <new>
#include <random>

void *operator new[](unsigned long size, char const *, int, unsigned int,
                     char const *, int) {
  return ::operator new(size);
}

void *operator new[](unsigned long size, unsigned long align, unsigned long,
                     char const *, int, unsigned int, char const *, int) {
  return ::operator new(size, std::align_val_t(align));
}

namespace indigox::utils {

  std::string ToUpper(const std::string &s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::toupper);
    return t;
  }

  std::string ToLower(const std::string &s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    return t;
  }

  std::string ToUpperFirst(const std::string &s) {
    std::string t = s;
    std::transform(t.begin(), t.begin() + 1, t.begin(), ::toupper);
    std::transform(t.begin() + 1, t.end(), t.begin() + 1, ::tolower);
    return t;
  }

  std::string GetRandomString(size_t length, size_t seed) {
    static std::string chrs =
        "qwertyuiopasdfghjklzxcvbnmZAQXSWCDEVFRBGTNHYMJUKILOP";
    static std::mt19937 rg{std::random_device{}()};
    static std::uniform_int_distribution<size_t> pick(0, chrs.size() - 1);
    if (seed != 0)
      rg.seed(seed);

    std::string s;
    s.reserve(length + 4); // no need to copy when adding extensions

    while (length--)
      s += chrs[pick(rg)];
    return s;
  }

} // namespace indigox::utils
