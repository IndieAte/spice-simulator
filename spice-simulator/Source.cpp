#include <iostream>
#include "Component.h"

int main() {

	std::vector<Component*> components;

	Resistor r("R1", 100, 0, 1);

	components.push_back(&r);

	std::cout << components[0]->getConductance(0, 1) << std::endl;

}