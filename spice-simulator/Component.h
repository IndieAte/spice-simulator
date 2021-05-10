#pragma once

struct Node;

#include "Node.h"
#include <string>

// Component Structure
// Base structure for components to be inhertited from
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
