#pragma once

#include "ConductanceMatrix.h"
#include <iostream>

std::vector<Eigen::Vector3d> runACAnalysis(int outputNode, double startFreq, double stopFreq, int pointsPerDecade,
	std::vector<Component*> components, int highest_node);