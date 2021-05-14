#pragma once

#include <iostream>
#include <complex>
#include <cmath>
#include "Eigen/Dense"
#include "Component.h"

std::vector<Eigen::Vector3d> runACAnalysis(int outputNode, double startFreq, double stopFreq, int pointsPerDecade,
	std::vector<Component*> components, int highest_node);