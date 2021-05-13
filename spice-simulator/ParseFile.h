#pragma once

#include <iostream>
#include <fstream>
#include "Component.h"
#include "SimulationCommand.h"

std::vector<Component*> decode_file(std::ifstream& infile, int& n, Command*& command);
double radians_to_degrees(double d);