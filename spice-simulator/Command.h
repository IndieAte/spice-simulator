#pragma once

#include <string>
#include <vector>

struct Command {
  Command(std::string p_type) : type { p_type } {}
  
  std::string type;

  virtual std::vector<double> getValues() {
    std::vector<double> v;
    return v;
  };
};

struct ACCommand : Command {
  ACCommand(std::string p_type, double p_sweep, double p_number_of_points,
    double p_start_frequency, double p_stop_frequency) : Command { p_type }, 
    sweep { p_sweep }, number_of_points { p_number_of_points},
    start_frequency { p_start_frequency }, stop_frequency { p_stop_frequency } {}

  double sweep;
  double number_of_points;
  double start_frequency;
  double stop_frequency;

  std::vector<double> getValues() {
    std::vector<double> v;
    v.push_back(sweep);
    v.push_back(number_of_points);
    v.push_back(start_frequency);
    v.push_back(stop_frequency);
    return v;
  }
};

struct OPCommand : Command {
  OPCommand(std::string p_type) : Command { p_type } {}

  std::vector<double> getValues() {
    std::vector<double> v;
    return v;
  }
};