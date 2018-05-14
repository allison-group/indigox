#include <algorithm>
#include <random>

#include <indigox/utils/common.hpp>
#include <indigox/utils/numerics.hpp>

namespace indigox::utils {
  
  string_ ToUpper(const string_& s){
    string_ t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::toupper);
    return t;
  }
  
  string_ ToLower(const string_& s){
    string_ t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    return t;
  }
  
  string_ ToUpperFirst(const string_& s){
    string_ t = s;
    std::transform(t.begin(), t.begin() + 1, t.begin(), ::toupper);
    std::transform(t.begin() + 1, t.end(), t.begin() + 1, ::tolower);
    return t;
  }
  
  string_ GetRandomString(size_ length, size_ seed) {
    static string_ chrs = "qwertyuiopasdfghjklzxcvbnmZAQXSWCDEVFRBGTNHYMJUKILOP";
    static std::mt19937 rg{std::random_device{}()};
    static std::uniform_int_distribution<size_t> pick(0, chrs.size() - 1);
    if (seed != 0) rg.seed(seed);
    
    string_ s;
    s.reserve(length + 4);  // no need to copy when adding extensions
    
    while(length--)
      s += chrs[pick(rg)];
    return s;
  }
  
}  // namespace indigox::utils
