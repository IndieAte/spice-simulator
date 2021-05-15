#pragma once

#include <iostream>

#include "Component.h"
#include "Eigen/Dense"

Eigen::VectorXd runDCOpPoint(std::vector<Component*> comps, int nNodes);