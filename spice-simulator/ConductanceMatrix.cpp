#include "ConductanceMatrix.h"

using namespace Eigen;

void updateConductanceMatrix(MatrixXcd& conductanceMatrix, Component* component) {
	if (typeid(*component) == typeid(CurrentSource)) {
		return;
	} else {
		std::vector<int> nodes = component->getNodes();

		if (nodes.size() == 2) {
			if (nodes[0] != 0 && nodes[1 != 0]) {
				conductanceMatrix(nodes[0] - 1, nodes[1] - 1) -= component->getConductance(nodes[0], nodes[1]);
				conductanceMatrix(nodes[1] - 1, nodes[0] - 1) -= component->getConductance(nodes[1], nodes[0]);
			}

			if (nodes[0] != 0) {
				conductanceMatrix(nodes[0] - 1, nodes[0] - 1) += component->getConductance(nodes[0], nodes[1]);
			}

			if (nodes[1] != 0) {
				conductanceMatrix(nodes[1] - 1, nodes[1] - 1) += component->getConductance(nodes[1], nodes[0]);
			}
		}
	}
}

MatrixXcd getConductanceMatrix(std::vector<Component*> components, int numNodes) {
	MatrixXcd conductanceMatrix = MatrixXcd::Zero(numNodes, numNodes);

	for (int i = 0; i < components.size(); i++) {
		updateConductanceMatrix(conductanceMatrix, components[i]);
	}

	return conductanceMatrix;
}

VectorXcd getCurrentVector(std::vector<Component*> components, int numNodes) {
	VectorXcd currentVector(numNodes);

	return currentVector;
}