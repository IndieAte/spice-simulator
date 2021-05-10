#include "ComplexMatrix.h"
#include "Component.h"

using namespace Eigen;

ComplexMatrix getConductanceMatrix(std::vector<Component> components, int numNodes) {
	ComplexMatrix conductionMatrix(numNodes, numNodes);

	return conductionMatrix;
}

ComplexVector getCurrentVector(std::vector<Component> components, int numNodes) {
	ComplexVector currentVector(numNodes);

	return currentVector;
}