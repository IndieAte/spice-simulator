#pragma once
#include "Node.h"
#include <vector>
#include <string>

struct Component {
  std::string name;
  std::vector<Node*> nodes;
};

struct CurrentSource : Component {
  double current;
};

struct Resistor : Component {
  double resistance;
};