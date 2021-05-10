#include <iostream>
#include "ConductanceMatrix.h"
#include "Component.h"

int main() {
	std::vector<Component*> components;

	components.push_back(new CurrentSource("I1", 0.002, 2, 3));
	components.push_back(new CurrentSource("I2", 0.001, 4, 2));


	components.push_back(new Resistor("R1", 1000, 0, 1));
	components.push_back(new Resistor("R2", 2000, 1, 2));
	components.push_back(new Resistor("R3", 1000, 3, 0));
	components.push_back(new Resistor("R4", 3000, 4, 3));

	int numNodes = 4;

	Eigen::VectorXcd solution = solveAtFrequency(components, numNodes, 0);

	for (int i = 0; i < numNodes; i++) {
		std::cout << "Node " << i + 1 << ": ";
		std::cout << real(solution(i)) << "V" << std::endl;
	}

	for (int i = 0; i < components.size(); i++) {
		delete(components[i]);
	}
}