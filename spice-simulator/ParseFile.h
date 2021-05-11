#pragma once

#include <iostream>
#include <fstream>
#include "Component.h"

std::vector<Component*> decode_file(std::ifstream& infile, int& n);