#pragma once

#include <string>
#include <vector>

// Component Class
// Base class for component object, with basic constructor to set component name
// and the virtual prototypes for two functions, getNodes and getConductance
class Component {
public:
  Component(std::string p_name) : name{ p_name } {}

  // getNodes Function
  // Implemented by each derived class, returns a vector of the integer IDs of
  // each node connected to the component
  virtual std::vector<int> getNodes();

  // getConductance Function
  // Implemented by each derived class, returns the conductance between two
  // nodes, so long as the component connects those nodes
  virtual double getConductance(int p_node1, int p_node2);

  // getProperties Function
  // Implemented by each derived class, returns a vector of component properties
  // to describe it's behaviour
  virtual std::vector<double> getProperties();

protected:
  std::string name;
};

// CurrentSource Class
// Derived from Component, implements a DC current source
class CurrentSource : public Component {
public:
  CurrentSource(std::string p_name, double p_current, int p_nodeIn, int p_nodeOut) : 
    Component{ p_name }, current{ p_current }, 
    nodeIn{ p_nodeIn }, nodeOut{ p_nodeOut } {}

  std::vector<int> getNodes() override;
  double getConductance(int p_node1, int p_node2) override;
  std::vector<double> getProperties() override;

private:
  double current;
  int nodeIn, nodeOut;
};

// Resistor Class
// Derived from Component, implements a resistor
class Resistor : public Component{
public:
  Resistor(std::string p_name, double p_resistance, int p_node1, int p_node2) : 
    Component{ p_name }, resistance{ p_resistance }, node1{ p_node1 },
    node2{ p_node2 } {}

  std::vector<int> getNodes() override;
  double getConductance(int p_node1, int p_node2) override;
  std::vector<double> getProperties() override;

private:
  double resistance;
  int node1, node2;
};
