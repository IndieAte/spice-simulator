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

	Eigen::MatrixXcd conductanceMatrix = getConductanceMatrix(components, 2);
	Eigen::VectorXcd currentVector = getCurrentVector(components, 2);

	Eigen::VectorXcd voltageVector;

	voltageVector = conductanceMatrix.colPivHouseholderQr().solve(currentVector);

	std::cout << voltageVector << std::endl;
}