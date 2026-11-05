#pragma once

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace worm_potential {

// Dense site order: basis + n_basis * (x + Lx * (y + Ly * z)).
// Whitespace separates values; '#' starts a comment anywhere on a line.
inline std::vector<double> read_dense(const std::string& filename, std::size_t nsites) {
  std::ifstream input(filename);
  if (!input) throw std::runtime_error("Cannot open onsite_potential_file: " + filename);
  std::vector<double> values;
  values.reserve(nsites);
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    std::istringstream fields(line);
    std::string token;
    while (fields >> token) {
      std::size_t consumed = 0;
      double value;
      try { value = std::stod(token, &consumed); }
      catch (const std::exception&) {
        throw std::runtime_error("Invalid potential value at " + filename + ":" + std::to_string(line_number));
      }
      if (consumed != token.size() || !std::isfinite(value))
        throw std::runtime_error("Potential values must be finite real numbers at " + filename + ":" + std::to_string(line_number));
      if (values.size() == nsites)
        throw std::runtime_error("Too many values in onsite_potential_file (expected " + std::to_string(nsites) + ").");
      values.push_back(value);
    }
  }
  if (input.bad()) throw std::runtime_error("Failed reading onsite_potential_file: " + filename);
  if (values.size() != nsites)
    throw std::runtime_error("onsite_potential_file has " + std::to_string(values.size()) +
                             " values; expected " + std::to_string(nsites) + ".");
  return values;
}

} // namespace worm_potential
