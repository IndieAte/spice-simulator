#include "ConductanceMatrix.h"
#include "Component.h"
#include "ParseFile.h"

int main() {
	std::vector<Component*> components;

	components.push_back(new ACVoltageSource("V1", 1, 0, 1, 0));
	components.push_back(new Resistor("R1", 1000, 1, 2));
	components.push_back(new Resistor("R2", 1000, 2, 0));

	int highestNode = 2;

	Eigen::MatrixXcd GMat = getConductanceMatrix(components, highestNode);
	Eigen::VectorXcd IVec = getCurrentVector(components, highestNode);

	std::cout << "The conduction matrix is:" << std::endl;
	std::cout << GMat << std::endl;

	std::cout << "The current vector is:" << std::endl;
	std::cout << IVec << std::endl;

	Eigen::VectorXcd VVec = solveAtFrequency(components, highestNode, 100);

	std::cout << "The voltage vector is:" << std::endl;
	std::cout << VVec << std::endl;
}