#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include "Component.h"
#include "Command.h"
#include "Model.h"

std::vector<Component*> decode_file(std::ifstream& infile, int& n, Command*& command, std::vector<int>& ac_source_indexes);
double radians_to_degrees(double d);
bool is_number(std::string s, bool dec_check);