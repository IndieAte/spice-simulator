#pragma once

#include <complex>
#include <Eigen/Dense>
#include "Component.h"

Eigen::MatrixXcd getConductanceMatrix(std::vector<Component*> components, int numNodes);