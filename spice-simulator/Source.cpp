#include <iostream>
#include "Component.h"

int main() {

	Resistor r("R1", 100, 0, 1);

	std::cout << r.getConductance(0, 1) << std::endl;

}