#pragma once
#include "Node.h"
#include <vector>

struct Component {
  std::vector<Node*> nodes;
};

struct CurrentSource : Component {
  double current;
};

struct Resistor : Component {
  double resistance;
};