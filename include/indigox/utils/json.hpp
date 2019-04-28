#include <nlohmann/json.hpp>

#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifndef INDIGOX_UTILS_JSON_HPP
#define INDIGOX_UTILS_JSON_HPP

using json = nlohmann::json;

namespace indigox::utils {
  std::string JSONKeyChecker(std::vector<std::string> &required_keys,
                             std::vector<std::string> &optional_keys,
                             json &data) {
    std::set<std::string> req(required_keys.begin(), required_keys.end());
    std::set<std::string> opt(optional_keys.begin(), optional_keys.end());
    std::stringstream error_stream;

    bool added_additional = false;

    for (auto &[key, value] : data.items()) {
      auto req_pos = req.find(key);
      auto opt_pos = opt.find(key);
      if (req_pos != req.end())
        req.erase(req_pos);
      else if (opt_pos == opt.end()) {
        if (!added_additional) {
          error_stream << "Additional keys=";
          added_additional = true;
        }
        error_stream << key << ",";
      }
    }

    if (!req.empty()) {
      error_stream << "Missing required keys=";
      for (auto &key : req) error_stream << key << ",";
    }

    return error_stream.str();
  }
} // namespace indigox::utils

#endif /* INDIGOX_UTILS_JSON_HPP */
