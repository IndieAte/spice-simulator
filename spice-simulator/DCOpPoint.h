#pragma once

#include "Component.h"
#include "Eigen/Dense"

Eigen::VectorXcd runDCOpPoint(std::vector<Component*> comps, int nNodes);