#include "ConductanceMatrix.h"

using namespace Eigen;

// updateConductanceMatrix Function
// Updates the conductance matrix (passed by reference) with the effect of each
// component as this function is called for it
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

// updateCurrentVector Function
// Helper function that updates the current vector (passed by reference), which is only
// affected by sources and non-linear components
void updateCurrentVector(VectorXcd& currentVector, Component* component) {
	if (typeid(*component) == typeid(CurrentSource)) {
		// Input at 0, Output at 1
		std::vector<int> nodes = component->getNodes();
		double current = (component->getProperties())[0];

		if (nodes[0] != 0) {
			currentVector(nodes[0] - 1) -= current;
		}

		if (nodes[1] != 0) {
			currentVector(nodes[1] - 1) += current;
		}
	}
}

// getConductanceMatrix Function
// Returns a matrix of complex doubles representing the conductance matrix in the
// nodal analysis equation
MatrixXcd getConductanceMatrix(std::vector<Component*> components, int numNodes) {
	MatrixXcd conductanceMatrix = MatrixXcd::Zero(numNodes, numNodes);

	for (int i = 0; i < components.size(); i++) {
		updateConductanceMatrix(conductanceMatrix, components[i]);
	}

	return conductanceMatrix;
}

// getCurrentVector Function
// Returns a vector of complex doubles representing the current vector in the
// nodal analysis equation
VectorXcd getCurrentVector(std::vector<Component*> components, int numNodes) {
	VectorXcd currentVector = VectorXcd::Zero(numNodes);

	for (int i = 0; i < components.size(); i++) {
		updateCurrentVector(currentVector, components[i]);
	}

	return currentVector;
}

// solveAtFrequency Function
// Solves for the voltage vector given a vector of pointers to Components, the number
// of nodes in the circuit and the frequency to solve for
VectorXcd solveAtFrequency(std::vector<Component*> components, int numNodes, double frequency) {
	MatrixXcd conductanceMatrix = getConductanceMatrix(components, numNodes);
	VectorXcd currentVector = getCurrentVector(components, numNodes);

	VectorXcd solution(numNodes);

	solution = conductanceMatrix.colPivHouseholderQr().solve(currentVector);

	return solution;
}