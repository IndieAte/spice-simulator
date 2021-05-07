#pragma once
#include "Node.h"
#include <vector>
#include <string>

struct Component {
  std::string name;
};

struct CurrentSource : Component {
  double current;

  Node* in;
  Node* out;
};

struct Resistor : Component {
  double resistance;

  Node* node1;
  Node* node2;
};