#include "ComplexMatrix.h"

using namespace Eigen;

void updateConductionMatrix(MatrixXcd& conductionMatrix, Component* component);

MatrixXcd getConductanceMatrix(std::vector<Component*> components, int numNodes) {
	MatrixXcd conductionMatrix = MatrixXcd::Zero(numNodes, numNodes);

	for (int i = 0; i < components.size(); i++) {
		updateConductionMatrix(conductionMatrix, components[i]);
	}

	return conductionMatrix;
}

ComplexVector getCurrentVector(std::vector<Component*> components, int numNodes) {
	ComplexVector currentVector(numNodes);

	return currentVector;
}

void updateConductionMatrix(MatrixXcd& conductionMatrix, Component* component) {
	if (typeid(*component) == typeid(CurrentSource)) {
		return;
	} else {
		std::vector<int> nodes = component->getNodes();

		if (nodes.size() == 2) {
			if (nodes[0] != 0 && nodes[1 != 0]) {
				conductionMatrix(nodes[0] - 1, nodes[1] - 1) += component->getConductance(nodes[0], nodes[1]);
				conductionMatrix(nodes[1] - 1, nodes[0] - 1) += component->getConductance(nodes[1], nodes[0]);
			}

			if (nodes[0] != 0) {
				conductionMatrix(nodes[0] - 1, nodes[0] - 1) += component->getConductance(nodes[0], nodes[1]);
			}

			if (nodes[1] != 0) {
				conductionMatrix(nodes[1] - 1, nodes[1] - 1) += component->getConductance(nodes[1], nodes[0]);
			}
		}
	}
}