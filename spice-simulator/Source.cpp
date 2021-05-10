#include <iostream>
#include "ConductanceMatrix.h"
#include "Component.h"

int main() {
	std::vector<Component*> components;

	CurrentSource i1("I1", 1, 0, 1);
	Resistor r1("R1", 100, 1, 2);
	Resistor r2("R2", 200, 2, 0);

	components.push_back(&i1);
	components.push_back(&r1);
	components.push_back(&r2);

	std::cout << getConductanceMatrix(components, 2) << std::endl;
}