#include "DCOpPoint.h"

using namespace Eigen;

VectorXcd runDCOpPoint(std::vector<Component*> comps, int nNodes) {
	std::vector<int> cSIndexes, vSIndexes, lCIndexes, nlCIndexes;

	for (int i = 0; i < comps.size(); i++) {
		Component* c = comps[i];

		if (typeid(*c) == typeid(ACCurrentSource) || typeid(*c) == typeid(DCCurrentSource)) {
			cSIndexes.push_back(i);
		} else if (typeid(*c) == typeid(ACVoltageSource) || typeid(*c) == typeid(DCVoltageSource) ||
			typeid(*c) == typeid(Inductor)) {
			vSIndexes.push_back(i);
		} else if (typeid(*c) == typeid(Diode)) {
			nlCIndexes.push_back(i);
		} else {
			lCIndexes.push_back(i);
		}
	}

	VectorXcd prevSoln, currSoln;
	return currSoln;
}

VectorXcd iterate(std::vector<Component*> comps, std::vector<int> cSIndexes, std::vector<int> vSIndexes,
	std::vector<int> lCIndexes, std::vector<int> nlCIndexes, int nNodes) {
	
	MatrixXcd gMat = MatrixXcd::Zero(nNodes, nNodes);
	VectorXcd iVec = VectorXcd::Zero(nNodes);

	for (int i = 0; i < lCIndexes.size(); i++) {
		int j = lCIndexes[i];
		linearComponentHandler(comps[j], gMat, iVec);
	}

	for (int i = 0; i < nlCIndexes.size(); i++) {
		int j = nlCIndexes[i];
		nonlinearComponentHandler(comps[j], gMat, iVec);
	}

	for (int i = 0; i < cSIndexes.size(); i++) {
		int j = cSIndexes[i];
		currentSourceHandler(comps[j], iVec);
	}

	for (int i = 0; i < vSIndexes.size(); i++) {
		int j = vSIndexes[i];
		voltageSourceHandler(comps[j], gMat, iVec);
	}

	VectorXcd vVec = gMat.colPivHouseholderQr().solve(iVec);

	for (int i = 0; i < nlCIndexes.size(); i++) {
		int j = nlCIndexes[i];
		updateNonlinearComponent(comps[j], vVec);
	}
}

void linearComponentHandler(Component* comp, MatrixXcd& gMat, VectorXcd& iVec) {
	std::vector<int> nodes = comp->getNodes();

	int n0 = nodes[0];
	int n1 = nodes[1];

	if (n0 == n1) return;

	if (typeid(*comp) == typeid(Resistor)) {
		int n0i = n0 - 1;
		int n1i = n1 - 1;

		std::complex<double> g = comp->getConductance(n0, n1, 0);

		if (n0 != 0 && n1 != 0) {
			gMat(n0i, n1i) -= g;
			gMat(n1i, n0i) -= g;
		}

		if (n0 != 0) gMat(n0i, n0i) += g;
		if (n1 != 0) gMat(n1i, n1i) += g;
	}
}

void nonlinearComponentHandler(Component* comp, MatrixXcd& gMat, VectorXcd& iVec) {

}

void currentSourceHandler(Component* comp, VectorXcd& iVec) {
	if (typeid(*comp) == typeid(DCCurrentSource)) {
		std::vector<int> nodes = comp->getNodes();
		int nIn = nodes[0] - 1;
		int nOut = nodes[1] - 1;

		if (nIn == nOut) return;

		std::vector<double> ppts = comp->getProperties();
		double current = ppts[0];

		if (nIn != -1 && nOut != -1) {
			iVec(nIn) -= current;
			iVec(nOut) += current;
		} else if (nIn == -1) {
			iVec(nOut) += current;
		} else {
			iVec(nIn) -= current;
		}
	}
}

void voltageSourceHandler(Component* comp, MatrixXcd& gMat, VectorXcd& iVec) {

}

void updateNonlinearComponent(Component* comp, VectorXcd vVec) {
	
}