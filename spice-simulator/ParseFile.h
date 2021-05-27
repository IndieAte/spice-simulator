#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include "Component.h"
#include "Command.h"
#include "Model.h"

std::vector<Component*> decode_file(std::ifstream& infile, int& n, Command*& command);
double radians_to_degrees(double d);