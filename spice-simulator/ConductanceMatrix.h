#pragma once

#include <complex>
#include "Eigen/Dense"
#include "Component.h"

Eigen::VectorXcd solveAtFrequency(std::vector<Component*> components, int numNodes, double frequency);