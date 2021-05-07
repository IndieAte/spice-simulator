#pragma once

#include <iostream>
#include <vector>
#include "Component.h"

// Node Structure
// Created for each circuit node, contains identifying information
// and a vector of pointers to connected components
struct Node {
	std::string name;
	int number;

	std::vector<Component*> components;
};