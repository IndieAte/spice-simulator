#pragma once

#include <complex>
#include "Eigen/Dense"
#include "Component.h"

Eigen::MatrixXcd getConductanceMatrix(std::vector<Component*> components, int numNodes);
Eigen::VectorXcd getCurrentVector(std::vector<Component*> components, int numNodes);
Eigen::VectorXcd solveAtFrequency(std::vector<Component*> components, int numNodes, double frequency);